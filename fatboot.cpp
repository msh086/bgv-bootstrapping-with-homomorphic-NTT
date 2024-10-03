/* Copyright (C) 2012-2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>
#include <helib/FHE.h>
#include <helib/EncryptedArray.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>
#include <helib/sample.h>
#include <helib/ArgMap.h>
#include <helib/EvalMap.h>


bool enable_lin_test = true;

long bit_reverse(long x, long nbits) {
    long res = 0;
    for(long i = 0; i < nbits; i++)
        if(x & (1L << i))
            res += 1L << (nbits - 1 - i);
    return res + ((x >> nbits) << nbits); // keep the higher bits
}

std::vector<long> auto_po2_parition(long m, long p, bool forceRadix2) {
    helib::assertTrue(helib::is2power(m), "???");
    long logm = NTL::NumTwos(NTL::ZZ(m));
    long nmats, A, logd, logD, logg;
    long split_pt;
    if (p % 4 == 1) {
        A = helib::ord(p-1, 2);
        logd = std::max(logm - A, 0L);
        logD = logm - 2 - logd;
        logg = (logD + 1) >> 1;
        // the diagonal indices of the first IFFT matrix is divisible by g
        split_pt = logD - logg;
    } else {
        A = helib::ord(p+1, 2);
        logd = std::max(logm - A, 1L);
        logD = logm - 1 - logd;
        logg = (logD + 1) >> 1;
        if (forceRadix2)
            split_pt = logD - logg;
        else
            split_pt = logD - logg + 1;
    }
    nmats = logm - 1 - logd + 1;
    // [0, split_pt] + [split_pt + 1, nmats]
    return {0, split_pt + 1, nmats};
}

// #include "gtest/gtest.h"
// #include "test_common.h"
struct Parameters
{
    const long p; // plaintext base
    const long r; // exponent
    // p^r is the plaintext-space modulus
    const long c;        // number of columns in the key-switching matrices
    const long bits;     // # of bits in the modulus chain
    const long skHwt;    // Hamming weight of recryption secret key
    const long nthreads; // number of threads
    const long seed;     // random number seed
    int useCache;  // 0: zzX cache, 1: DCRT cache
    const int force_bsgs;
    const int force_hoist;
    // int disable_intFactor // fhe_disable_intFactor (disabled by Victor)
    const int chen_han;
    const bool debug; // generate debugging output
    const int scale;  // scale parameter
    NTL::Vec<long> global_gens;
    NTL::Vec<long> global_ords;
    const NTL::Vec<long> global_mvec;
    const int c_m; // = 100;
    const long outer_rep;
    long inner_rep;

    Parameters(long p,
               long r,
               long c,
               long bits,
               long skHwt,
               long nthreads,
               long seed,
               long useCache,
               int c_m,
               int force_bsgs,
               int force_hoist,
               int chen_han,
               bool debug,
               int scale,
               const std::vector<long> &global_gens,
               const std::vector<long> &global_ords,
               const std::vector<long> &global_mvec,
               long outer_rep,
               long inner_rep) : p(p),
                                 r(r),
                                 c(c),
                                 bits(bits),
                                 skHwt(skHwt),
                                 nthreads(nthreads),
                                 seed(seed),
                                 useCache(useCache),
                                 force_bsgs(force_bsgs),
                                 force_hoist(force_hoist),
                                 chen_han(chen_han),
                                 debug(debug),
                                 scale(scale),
                                 global_gens(helib::convert<NTL::Vec<long>>(global_gens)),
                                 global_ords(helib::convert<NTL::Vec<long>>(global_ords)),
                                 global_mvec(helib::convert<NTL::Vec<long>>(global_mvec)),
                                 c_m(c_m),
                                 outer_rep(outer_rep),
                                 inner_rep(inner_rep)
    {
        if (global_gens.empty() || global_ords.empty() || global_mvec.empty())
            throw helib::LogicError("gens, ords, and mvec must be non-empty");
    };

    friend std::ostream &operator<<(std::ostream &os, const Parameters &params)
    {
        return os << "{"
                  << "p" << params.p << ","
                  << "r" << params.r << ","
                  << "c" << params.c << ","
                  << "bits" << params.bits << ","
                  << "skHwt" << params.skHwt << ","
                  << "nthreads" << params.nthreads << ","
                  << "seed" << params.seed << ","
                  << "useCache" << params.useCache << ","
                  << "force_bsgs" << params.force_bsgs << ","
                  << "force_hoist" << params.force_hoist << ","
                  << "chen_han" << params.chen_han << ","
                  << "debug" << params.debug << ","
                  << "scale" << params.scale << ","
                  << "global_gens" << params.global_gens << ","
                  << "global_ords" << params.global_ords << ","
                  << "global_mvec" << params.global_mvec << ","
                  << "c_m" << params.c_m << "}";
    }
};

class GTestFatboot
{
    void preContextSetup()
    {

        if (!noPrint)
            helib::fhe_stats = true;

        if (!noPrint)
        {
            std::cout << "*** GTestFatboot";
            if (helib::isDryRun())
                std::cout << " (dry run)";
            std::cout << ": p=" << p << ", r=" << r << ", bits=" << bits
                      << ", t=" << skHwt << ", c=" << c << ", m=" << m
                      << " mvec=" << mvec << ", gens=" << helib::vecToStr(gens)
                      << ", ords=" << helib::vecToStr(ords) << std::endl;
            std::cout << "Computing key-independent tables..." << std::flush;
        }
        helib::setTimersOn();
        helib::setDryRun(
            false); // Need to get a "real context" to test bootstrapping
        time = -NTL::GetTime();
    }

    static void setGlobals(int force_bsgs, int force_hoist, int chen_han)
    {
        helib::fhe_test_force_bsgs = force_bsgs;
        helib::fhe_test_force_hoist = force_hoist;
        helib::fhe_force_chen_han = chen_han;
    };

    void cleanupBootstrappingGlobals()
    {
        helib::fhe_test_force_bsgs = old_fhe_test_force_bsgs;
        helib::fhe_test_force_hoist = old_fhe_test_force_hoist;
        helib::fhe_force_chen_han = old_fhe_force_chen_han;
    }

    static void setSeedIfNeeded(const long seed)
    {
        if (seed)
            SetSeed(NTL::ZZ(seed));
    };

    static void checkPM(const long p, const long m)
    {
        helib::assertTrue(NTL::GCD(p, m) == 1, "GCD(p, m) == 1");
    }

    const int old_fhe_test_force_bsgs;
    const int old_fhe_test_force_hoist;
    const int old_fhe_force_chen_han;
    const long p;
    const long r;
    const long c;
    const long bits;
    const long skHwt;
    const long nthreads;
    const long seed;
    const long useCache;
    const int force_bsgs;
    const int force_hoist;
    const int chen_han;
    const bool debug;
    const int scale;
    const NTL::Vec<long> mvec;
    const std::vector<long> gens;
    const std::vector<long> ords;
    const int c_m; // = 100;
    const long outer_rep;
    const long inner_rep;
    const bool noPrint;
    const bool dry;
    // new params
    const long encapSkHwt;
    const long btsC;
    const long t;
    const double btsScale;
    const bool newBtsFlag;
    const bool s2cFirstFlag;
    const bool isThick;
    // po2 linear transform
    const bool forceRadix2;
    const std::vector<long> partition;
    const bool dummy;

    const long m, phim;
    double time;
    helib::Context context;

public:
    GTestFatboot(Parameters param,
                 bool noPrint,
                 bool dry,
                 long encapSkHwt,
                 long btsC,
                 long t,
                 long btsScale,
                 long newBtsFlag,
                 bool newKSFlag,
                 bool s2cFirstFlag,
                 bool isThick,
                 bool forceRadix2,
                 std::vector<long> partition,
                 bool dummy) : 
                                 old_fhe_test_force_bsgs(helib::fhe_test_force_bsgs),
                                 old_fhe_test_force_hoist(helib::fhe_test_force_hoist),
                                 old_fhe_force_chen_han(helib::fhe_force_chen_han),
                                 p((setGlobals(param.force_bsgs,
                                               param.force_hoist,
                                               param.chen_han),
                                    param.p)),
                                 r(param.r),
                                 c(param.c),
                                 bits(param.bits),
                                 skHwt(param.skHwt),
                                 nthreads((NTL::SetNumThreads(param.nthreads), param.nthreads)),
                                 seed((setSeedIfNeeded(param.seed), param.seed)),
                                 useCache(param.useCache),
                                 force_bsgs(param.force_bsgs),
                                 force_hoist(param.force_hoist),
                                 chen_han(param.chen_han),
                                 debug(param.debug),
                                 scale(param.scale),
                                 mvec(param.global_mvec),
                                 gens(helib::convert<std::vector<long>>(param.global_gens)),
                                 ords(helib::convert<std::vector<long>>(param.global_ords)),
                                 c_m(param.c_m),
                                 outer_rep(param.outer_rep),
                                 inner_rep(param.inner_rep),
                                 noPrint(noPrint),
                                 dry(dry),
                                 encapSkHwt(encapSkHwt),
                                 btsC(btsC),
                                 t(t),
                                 btsScale(btsScale),
                                 newBtsFlag(newBtsFlag),
                                 s2cFirstFlag(s2cFirstFlag),
                                 isThick(isThick),
                                 forceRadix2(forceRadix2),
                                 partition(partition),
                                 dummy(dummy),
                                 m(helib::computeProd(mvec)),
                                 phim((checkPM(p, m), helib::phi_N(m))),
                                 time(0),
                                 context((preContextSetup(),
                                          helib::ContextBuilder<helib::BGV>()
                                              .m(m)
                                              .p(p)
                                              .r(r)
                                              .gens(gens)
                                              .ords(ords)
                                              .scale(scale ? scale : 10 /*default is 10.*/)
                                              // new params
                                              .encapSkHwt(encapSkHwt)
                                              .btsScale(btsScale)
                                              .btsC(btsC)
                                              .t(t)
                                              .newKS(newKSFlag)
                                              .newBts(newBtsFlag)
                                              .s2cFirst(s2cFirstFlag)
                                              // remove the hack
                                              .bits(bits)
                                              .c(c)
                                              .bootstrappable(true) // XXX: debug
                                              .skHwt(skHwt)
                                              .mvec(mvec)
                                              .buildCache(useCache)
                                              .setThick(isThick)
                                              .forceRadix2(forceRadix2)
                                              .FFT_partition(partition)
                                              .dummy(dummy) // XXX: debug
                                              //    .buildModChain(false)
                                              .build()))
    {
    }

    // NOTE: generate parameter for power-of-2 bts
    static Parameters genPo2Param(long logm, long logd, long ptype, long r, bool force_chen_han) {
        helib::assertTrue(logd < logm, "logd too large");
        long p;
        std::vector<long> gens, ords;
        if(ptype == 1) { // p == 1 mod 4
            helib::assertTrue(logm - logd - 2 >= 1, "logd too large for p=1 mod 4");
            gens.push_back(-1);
            gens.push_back(5);
            ords.push_back(2);
            ords.push_back(1L << (logm - logd - 2));
            if(logd == 0) { // p = 2^{A}*k+1 for A >= m
                long step = 1L << logm;
                p = step + 1;
                while(!NTL::ProbPrime(p))
                    p += step;
            } else { // p = 2^{A}*k+1 for A = logm-logd and odd k
                long step = 1L << (logm-logd+1);
                p = step/2 + 1;
                while(!NTL::ProbPrime(p))
                    p += step;
            }
        } else { // p == 3 mod 4
            helib::assertTrue(logd > 0, "logd = 0 for p=3 mod 4");
            gens.push_back(5);
            ords.push_back(1L << (logm - logd - 1));
            if(logd == 1) { // p = 2^{A}*k-1 for A>=logm-1
                long step = 1L << (logm - 1);
                p = step - 1;
                while(!NTL::ProbPrime(p))
                    p += step;
            } else { // p = 2^{A}*k-1 for A = logm-logd and odd k
                long step = 1L << (logm - logd + 1);
                p = step/2 - 1;
                while(!NTL::ProbPrime(p))
                    p += step;
            }
        }
        return Parameters(
                // p, r, c, bits, h
                p, r, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                gens, ords, {1L << logm}, 1, 5);
    }

    void TearDown()
    {
        if (noPrint)
        {
            helib::printAllTimers();
        }
        if (helib::fhe_stats)
            helib::print_stats(std::cout);

        cleanupBootstrappingGlobals();
        helib::cleanupDebugGlobals();
    }

    void run()
    {
        // context.buildModChain(bits,
        //                       c,
        //                       /*willBeBootstrappable=*/true,
        //                       /*t=*/skHwt);

        if (!noPrint)
        {
            std::cout << "security=" << context.securityLevel() << std::endl;
            std::cout << "# small primes = " << context.getSmallPrimes().card()
                      << std::endl;
            std::cout << "# ctxt primes = " << context.getCtxtPrimes().card()
                      << std::endl;
            std::cout << "# bits in ctxt primes = "
                      << long(context.logOfProduct(context.getCtxtPrimes()) / log(2.0) +
                              0.5)
                      << std::endl;
            std::cout << "# special primes = " << context.getSpecialPrimes().card()
                      << std::endl;
            std::cout << "# bits in special primes = "
                      << long(context.logOfProduct(context.getSpecialPrimes()) /
                                  log(2.0) +
                              0.5)
                      << std::endl;
            std::cout << "scale=" << context.getScale() << std::endl;
        }

        // context.enableBootStrapping(mvec, useCache);
        time += NTL::GetTime();

        if (!noPrint)
        {
            std::cout << " done in " << time << " seconds" << std::endl;
            std::cout << "  e=" << context.getRcData().e
                      << ", e'=" << context.getRcData().ePrime
                      << ", t=" << context.getRcData().skHwt << "\n"
                      << "  ";
            context.printout();
        }
        helib::setDryRun(
            dry); // Now we can set the dry-run flag if desired

        long p2r = context.getAlMod().getPPowR();
        helib::BootBench meanBenchmarker;

        for (long numkey = 0; numkey < outer_rep; numkey++)
        { // test with 3 keys
            time = -NTL::GetTime();
            if (!noPrint)
                std::cout << "Generating keys, " << std::flush;
            helib::SecKey secretKey(context);
            helib::PubKey &publicKey = secretKey;
            secretKey.GenSecKey(); // A +-1/0 secret key
            helib::addSome1DMatrices(
                secretKey); // compute key-switching matrices that we need
            helib::addFrbMatrices(secretKey);
            if (!noPrint)
                std::cout << "computing key-dependent tables..." << std::flush;
            // XXX: debug
            secretKey.genRecryptData();
            time += NTL::GetTime();
            if (!noPrint)
                std::cout << " done in " << time << " seconds\n";

            NTL::ZZX ptxt_poly1, ptxt_poly;
            std::vector<NTL::ZZX> val1;
            if (isThick)
            {
                NTL::zz_p::init(p2r);
                NTL::zz_pX poly_p = NTL::random_zz_pX(context.getPhiM());
                helib::zzX poly_p1 = helib::balanced_zzX(poly_p);
                ptxt_poly = helib::convert<NTL::ZZX>(poly_p1);
                helib::PolyRed(ptxt_poly1, ptxt_poly, p2r, true);
            }
            else
            {
                long d = context.getOrdP();
                long phim = context.getPhiM();
                long nslots = phim / d;

                NTL::zz_p::init(p2r);
                NTL::Vec<NTL::zz_p> val0(NTL::INIT_SIZE, nslots);

                bool debug_ordered = false;
                if (debug_ordered)
                    for (long i = 0; i < nslots; i++)
                        val0[i] = i + 1;
                else
                    for (auto &x : val0)
                        random(x);

                val1.resize(nslots);
                for (long i = 0; i < nslots; i++)
                    val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
            }
            // this is the format produced by decryption

            helib::setupDebugGlobals(&secretKey, context.shareEA());

            NTL::ZZX poly2;
            std::vector<NTL::ZZX> val2;
            helib::Ctxt c1(publicKey);

            if (isThick)
                secretKey.Encrypt(c1, ptxt_poly, p2r);
            else
                context.shareEA()->encrypt(c1, publicKey, val1);

#if 0
            // XXX: linear test
            bool is_po2 = helib::is2power(m);
            if (is_po2 && enable_lin_test)
            {
                auto ea = context.shareEA();
                long d = context.getOrdP();
                long nslots = context.getNSlots();
             
                std::cout << "start lin test\n";
                // test linear transform
                if (is_po2 && m <= 32)
                    for(auto factor : ea->getAlMod().getFactorsOverZZ())
                        std::cout << factor << "\n";
                // bool forceRadix2 = true;

                long nmats = 0;
                if(is_po2) {
                    std::cout << "power of 2 linear transform!\n";
                    // for thick bts, nmats = number of FFT layers + 1 (the encoding layer)
                    // although the encoding layer may be the identity map in thin bts,
                    // we still count it in for better compatibility
                    nmats = NTL::NumTwos(NTL::ZZ(context.getNSlots())) 
                        + 1; // log2(nslots) + 1
                    std::cout << "nmats = " << nmats << "\n";
                }
                // std::vector<long> partition{0, 4, nmats};
                helib::Ctxt linctxt = c1;
                if(isThick)
                    context.shareEA()->decrypt(linctxt, secretKey, val1);
                // thick: the coeff of lintest_poly will be moved into slots
                // thin:  the slots will be moved into coeffs
                NTL::ZZX lintest_poly;
                if(is_po2) {
                    secretKey.Decrypt(lintest_poly, linctxt);
                    lintest_poly.SetLength(context.getPhiM());
                    if(m <= 32) {
                        std::cout << "\n";
                        for(long i = 0; i < context.getPhiM(); i++)
                            std::cout << "poly[" << i << "] = " << NTL::coeff(lintest_poly, i) << '\n';
                        std::cout << "\n";
                        for(long i = 0; i < context.getNSlots(); i++)
                            for (long j = 0; j < ea->getDegree(); j++)
                                std::cout << "slot[" << i << "][" << j << "] = " << NTL::coeff(val1[i], j) << "\n";
                    }
                }
                // first map
                if(!is_po2){
                    // NOT UPDATED
                    auto firstMap = std::make_shared<helib::EvalMap>(*ea, false, mvec, true, useCache);
                    firstMap->apply(linctxt);
                } else {
                    // thick: IFFT (inv = 0)
                    // thin:  FFT (inv = 1)
                    if (isThick) {
                        auto firstMap = std::make_shared<helib::Po2IFFT>(*ea, false,
                            partition, forceRadix2, useCache);
                        firstMap->apply(linctxt);
                    } else {
                        auto firstMap = std::make_shared<helib::ThinPo2IFFT>(*ea, true,
                            partition, forceRadix2, useCache);
                        firstMap->apply(linctxt);
                    }
                }
                // thick: get the normal basis coefficients
                // thin:  get the polynomial coefficients
                context.getEA().restoreContext();
                std::vector<NTL::ZZX> vallin;
                NTL::vec_ZZ coeff_vec;
                NTL::vec_zz_p coeff_vec_p;
                context.getEA().decrypt(linctxt, secretKey, vallin);
                NTL::ZZX polylin;
                secretKey.Decrypt(polylin, linctxt);
                // for(long i = 0; i < context.getNSlots(); i++)
                //     for (long j = 0; j < ea->getDegree(); j++)
                //         std::cout << "slot[" << i << "] = " << NTL::coeff(vallin[i], j) << "\n";

                if(isThick) {
                    const NTL::mat_zz_p &p2n = 
                        context.getEA().getDerived(helib::PA_zz_p()).getNormalBasisMatrixInverse();
                    for(NTL::ZZX & ele : vallin) {
                        NTL::VectorCopy(coeff_vec, ele, d);
                        NTL::conv(coeff_vec_p, coeff_vec);
                        NTL::mul(coeff_vec_p, coeff_vec_p, p2n);
                        NTL::conv(coeff_vec, coeff_vec_p);
                        ele.rep = coeff_vec;
                    }
                }
                // thick: vallin vs lintest_poly
                // thin:  val1   vs polylin
                if(is_po2){
                    NTL::ZZX &target_poly = isThick ? lintest_poly : polylin;
                    std::vector<NTL::ZZX> &target_val = isThick ? vallin : val1;

                    if(m <= 32) { // debug print
                        std::cout << "\n";
                        if (isThick)
                            for(long i = 0; i < vallin.size(); i++) {
                                for(long j = 0; j < d; j++)
                                    std::cout << "slot[" << i << "][" << j << "] = " << 
                                        NTL::coeff(vallin[i], j) << "\n";
                            }
                        else
                            for(long i = 0; i < phim; i++)
                                std::cout << "poly[" << i << "] = " << NTL::coeff(polylin, i) << "\n";
                    }
                    // check bit reversal
                    if(p % 4 == 3){
                        if(!forceRadix2){
                            // unit = slot, all bits are reversed
                            long brbits = NTL::NumTwos(NTL::ZZ(nslots));
                            for(long i = 0; i < nslots; i++) {
                                for(long j = 0; j < d; j++) {
                                    helib::assertEq(NTL::coeff(target_val[bit_reverse(i, brbits)], j), 
                                        target_poly[i*d+j], "first Po2 map for p=3 mod 4 failed\n");
                                    if (!isThick && j > 0) {
                                        NTL::zz_p::init(p2r);
                                        NTL::zz_p tmp_rand;
                                        NTL::random(tmp_rand);
                                        NTL::conv(target_poly[i*d+j], tmp_rand);
                                    }
                                }
                            }
                            std::cout << "Bruun Po2 Step1 map for p=3 mod 4 ok\n";
                        } else {
                            // unit = half slot, all bits are reversed
                            long brbits = NTL::NumTwos(NTL::ZZ(nslots)) + 1;
                            for(long i = 0; i < nslots; i++) {
                                for(long j = 0; j < d; j++) {
                                    long dst_half_slot = bit_reverse(i*2 + j / (d/2), brbits);
                                    helib::assertEq(NTL::coeff(target_val[i], j),
                                        target_poly[dst_half_slot * d/2 + (j % (d/2))],
                                        "first forceRadix2 Po2 map for p=3 mod 4 failed\n");
                                    if (!isThick && j > 0) {
                                        NTL::zz_p::init(p2r);
                                        NTL::zz_p tmp_rand;
                                        NTL::random(tmp_rand);
                                        NTL::conv(target_poly[dst_half_slot * d/2 + (j % (d/2))], tmp_rand);
                                    }
                                }
                            }
                            std::cout << "first forceRadix2 Po2 map for p=3 mod 4 ok\n";
                            // XXX: removed
                            // NOTE: the coeffs are repacked into normal basis in this case
                            //  thus, we undo the P->N basis switch
                            // context.getEA().decrypt(linctxt, secretKey, vallin);
                        }
                    }
                    else { // p = 1 mod 4
                        // unit = slot, all but the highest bit are reversed
                        long brbits = NTL::NumTwos(NTL::ZZ(nslots)) - 1;
                        for(long i = 0; i < nslots; i++) {
                            for(long j = 0; j < d; j++) {
                                helib::assertEq(NTL::coeff(target_val[bit_reverse(i, brbits)], j), 
                                    target_poly[i*d+j], "first Po2 map for p=1 mod 4 failed\n");
                                if (!isThick && j > 0) {
                                    NTL::zz_p::init(p2r);
                                    NTL::zz_p tmp_rand;
                                    NTL::random(tmp_rand);
                                    NTL::conv(target_poly[i*d+j], tmp_rand);
                                }
                            }
                        }
                        std::cout << "Radix-2 Po2 Step1 map for p=1 mod 4 ok\n";
                    }
                }
#if 0
                // real unpacking from recryption.cpp
                std::vector<helib::Ctxt> unpacked(d, helib::Ctxt(helib::ZeroCtxtLike, linctxt));
                { // explicit scope to force all temporaries to be released
                    // prepare unpackSlotEncoding
                    const NTL::Mat<NTL::zz_p>& CBi =
                        ea->getDerived(helib::PA_zz_p()).getNormalBasisMatrixInverse();

                    std::vector<NTL::ZZX> LM;
                    LM.resize(d);
                    for (long i = 0; i < d; i++) // prepare the linear polynomial
                        LM[i] = rep(CBi[i][0]);

                    std::vector<NTL::ZZX> C;
                    ea->buildLinPolyCoeffs(C, LM); // "build" the linear polynomial

                    std::vector<NTL::ZZX> unpackSlotEncoding;
                    unpackSlotEncoding.resize(d); // encode the coefficients
                    // the coefficients above acts on a single slot, we need to make them SIMD
                    for (long j = 0; j < d; j++) {
                        std::vector<NTL::ZZX> v(nslots);
                        for (long k = 0; k < nslots; k++)
                        v[k] = C[j];
                        ea->encode(unpackSlotEncoding[j], v);
                    }   
                    // start unpacking
                    std::vector<std::shared_ptr<helib::DoubleCRT>> coeff_vector;
                    std::vector<double> coeff_vector_sz;
                    coeff_vector.resize(d);
                    coeff_vector_sz.resize(d);

                    for (long i = 0; i < d; i++) {
                    coeff_vector[i] = std::make_shared<helib::DoubleCRT>(unpackSlotEncoding[i],
                                                                    linctxt.getContext(),
                                                                    linctxt.getPrimeSet());
                    coeff_vector_sz[i] = NTL::conv<double>(
                        helib::embeddingLargestCoeff(unpackSlotEncoding[i],
                                                linctxt.getContext().getZMStar()));
                    }
                    std::vector<helib::Ctxt> frob(d, helib::Ctxt(helib::ZeroCtxtLike, linctxt));

                    NTL_EXEC_RANGE(d, first, last)
                    // FIXME: implement using hoisting!
                    for (long j = first; j < last; j++) { // process jth Frobenius
                    frob[j] = linctxt;
                    frob[j].frobeniusAutomorph(j);
                    frob[j].cleanUp();
                    // FIXME: not clear if we should call cleanUp here
                    }
                    NTL_EXEC_RANGE_END

                    helib::Ctxt tmp1(helib::ZeroCtxtLike, linctxt);
                    for (long i = 0; i < d; i++) {
                        for (long j = 0; j < d; j++) {
                            tmp1 = frob[j];
                            tmp1.multByConstant(*coeff_vector[helib::mcMod(i + j, d)],
                                                coeff_vector_sz[helib::mcMod(i + j, d)]);
                            unpacked[i] += tmp1;
                        }
                    }
                }
                // end of real unpacking
                for(long i = 0; i < d; i++) {
                    std::vector<NTL::ZZX> tmp_vec;
                    context.shareEA()->decrypt(unpacked[i], secretKey, tmp_vec);
                    for(long j = 0; j < tmp_vec.size(); j++) {
                        if(NTL::coeff(tmp_vec[j], 0) != NTL::coeff(vallin[j], i)){
                            std::cout << NTL::coeff(tmp_vec[j], 0) << " ### "
                                << NTL::coeff(vallin[j], i) << '\n';
                            exit(1);
                        }
                    }
                }
                std::cout << "unpack ok\n";
#endif
                if (isThick)
                    context.shareEA()->encrypt(linctxt, vallin);
                else
                    secretKey.Encrypt(linctxt, polylin, p2r);
                // second map
                if(!is_po2){
                    // NOT UPDATED
                    auto secondMap = std::make_shared<helib::EvalMap>(*ea,
                                                            false,
                                                            mvec,
                                                            false,
                                                            useCache);
                    secondMap->apply(linctxt);
                } else {
                    // thick: FFT (inv = 1)
                    // thin:  IFFT (inv = 0)
                    if (isThick) {
                        auto secondMap = std::make_shared<helib::Po2IFFT>(*ea, true, 
                            partition, forceRadix2, useCache);
                        secondMap->apply(linctxt);
                    } else {
                        auto secondMap = std::make_shared<helib::ThinPo2IFFT>(*ea, false,
                            partition, forceRadix2, useCache);
                        secondMap->apply(linctxt);
                    }
                }
                context.getEA().decrypt(linctxt, secretKey, vallin);
                helib::assertEq(vallin, val1, "lintest for native BTS failed");
                std::cout << "lintest ok\n";
            }
            // end linear test
            exit(0);
#endif
            NTL::SetNumThreads(1);
            helib::resetAllTimers();
            for (long num = 0; num < inner_rep; num++)
            { // multiple tests with same key
                if (isThick)
                {
                    helib::BootBench benchmarker;
                    if (t >= 0 && newBtsFlag)
                        benchmarker = publicKey.reCryptRefine(c1);
                    else // XXX: s2cFirst is not implemented for PubKey::reCrypt
                        benchmarker = publicKey.reCrypt(c1);
                    // print the current results
                    printf("time for linear1 = %f, linear2 = %f, extract = %f, total = %f\n",
                            benchmarker.time_linear_1, benchmarker.time_linear_2,
                            benchmarker.time_extract, benchmarker.time_total);
                    printf("bits for linear1 = %f, linear2 = %f, extract = %f, min cap = %f, after cap = %f, after prod = %f\n",
                            benchmarker.bits_down_linear_1, benchmarker.bits_down_linear_2,
                            benchmarker.bits_down_extract, context.getMinCap(), benchmarker.bits_final, benchmarker.bits_after_inner_prod);
                    printf("cap / time = %f\n",
                           (benchmarker.bits_final - context.getMinCap()) / benchmarker.time_total);
                    // add to the mean benchmarker
                    meanBenchmarker += benchmarker;
                    secretKey.Decrypt(poly2, c1);
                    ptxt_poly1.SetLength(context.getPhiM());
                    poly2.SetLength(context.getPhiM());
                    for (long i = 0; i < context.getPhiM(); i++)
                    {
                        if (ptxt_poly1[i] != poly2[i])
                            std::cout << i << " is bad\n";
                    }
                    helib::assertEq(ptxt_poly1, poly2, "fat results not match");
                }
                else
                {
                    helib::BootBench benchmarker;
                    if (t >= 0 && newBtsFlag)
                        benchmarker = publicKey.thinReCryptRefine(c1);
                    else
                        benchmarker = publicKey.thinReCrypt(c1);
                    printf("time for linear1 = %f, linear2 = %f, extract = %f, total = %f\n",
                            benchmarker.time_linear_1, benchmarker.time_linear_2,
                            benchmarker.time_extract, benchmarker.time_total);
                    printf("bits for linear1 = %f, linear2 = %f, extract = %f, min cap = %f, after cap = %f, after prod = %f\n",
                            benchmarker.bits_down_linear_1, benchmarker.bits_down_linear_2,
                            benchmarker.bits_down_extract, context.getMinCap(), benchmarker.bits_final, benchmarker.bits_after_inner_prod);
                    printf("cap / time = %f\n",
                           (benchmarker.bits_final - context.getMinCap() - benchmarker.bits_down_linear_1) / benchmarker.time_total);
                    meanBenchmarker += benchmarker;
                    context.shareEA()->decrypt(c1, secretKey, val2);
                    if (!dummy) { // no correctness check for dummy computation...
                        for (long i = 0; i < val1.size(); i++)
                        {
                            if (val1[i] != val2[i])
                                std::cout << i << " is bad\n";
                        }
                        helib::assertEq(val1, val2, "thin results not match");
                    } else
                        printf("warning: result not checked for dummy run\n");
                }
            }
        }
        meanBenchmarker.Mult(1.0 / (inner_rep * outer_rep));            
        printf("\n### summary over %ld test\n", inner_rep * outer_rep);        
        printf("time for linear1 = %f, linear2 = %f, extract = %f, total = %f\n",
                meanBenchmarker.time_linear_1, meanBenchmarker.time_linear_2,
                meanBenchmarker.time_extract, meanBenchmarker.time_total);
        printf("bits for linear1 = %f, linear2 = %f, extract = %f, min cap = %f, after cap = %f, after prod = %f\n",
                meanBenchmarker.bits_down_linear_1, meanBenchmarker.bits_down_linear_2,
                meanBenchmarker.bits_down_extract, context.getMinCap(), meanBenchmarker.bits_final, meanBenchmarker.bits_after_inner_prod);
        if (isThick)
            printf("cap / time = %f\n",
                           (meanBenchmarker.bits_final - context.getMinCap()) / meanBenchmarker.time_total);
        else
            printf("cap / time = %f\n",
                (meanBenchmarker.bits_final - context.getMinCap() - meanBenchmarker.bits_down_linear_1) / meanBenchmarker.time_total);

        printf("\n\n### bts finished, everything ok ###\n\n");
        if (!noPrint)
            helib::printAllTimers();
#if (defined(__unix__) || defined(__unix) || defined(unix))
        struct rusage rusage;
        getrusage(RUSAGE_SELF, &rusage);
        if (!noPrint)
            std::cout << "  rusage.ru_maxrss=" << rusage.ru_maxrss << std::endl;
#endif
        if (helib::fhe_stats)
            helib::print_stats(std::cout);
    }
};

int main(int argc, char *argv[])
{
    std::cout << std::unitbuf; 
    std::cerr << std::unitbuf;

    // NTL::zz_p::init(5);
    // NTL::mat_zz_p M0, M1;
    // M1.SetDims(2, 2);
    // M1 *= 2;
    // M0 *= 3;
    // auto M2 = M0 * M1;
    // if(NTL::IsZero(M0))
    //     printf("M0 empty\n");
    // M0 += M1;
    // std::cout << M0 << "\n";
    // return 0;


    helib::ArgMap amap;

    long i_arg = 0;
    amap.arg("i", i_arg, "index of the chosen parameter set");

    long h_arg = 0;
    amap.arg("h", h_arg, "hwt of encapsulated key");

    long t_arg = 0;
    amap.arg("t", t_arg, "parameter t used in new bts");

    long newBts_flag = 0;
    amap.arg("newbts", newBts_flag, "if new bts is used");

    long newKS_flag = 1;
    amap.arg("newks", newKS_flag, "if new ks is used");

    long thick_flag = 0;
    amap.arg("thick", thick_flag, "if thick bts is used");

    long repeat_arg = 5;
    amap.arg("repeat", repeat_arg, "number of tests");

    long forceRadix2_arg = 0;
    amap.arg("rad2", forceRadix2_arg, "force to use radix-2");

    long dummy_arg = 0;
    // amap.arg("dummy", dummy_arg, "dummy linear transform to simulate GV23");
    long baseline_arg = 0;
    amap.arg("baseline", baseline_arg, "to test baseline");

    long cache_arg = 1;
    amap.arg("cache", cache_arg, "use cache or not");

    long s2cFirst_flag = 0;
    amap.arg("s2cFirst", s2cFirst_flag, "if SlotToCoeff is the first step in general bootstrapping");

    amap.parse(argc, argv);
    // some candidate primes:
    // 2^14 - 3
    // 2^13 - 1
    // p,r,c, bits,skHwt,nthreads,seed,useCache,c_m,force_bsgs,force_hoist,chen_han,debug,scale,global_gens,global_ords,global_mvec,outer_rep, inner_rep
    bool force_chen_han = false;
    // clang-format off
    Parameters toy_params[] = {
        // NOTE: toy param sets, only for debug
        Parameters(
                // p, r, c, bits, h
                    17, 4, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {262, 146, 181}, {4, 6, 7}, {5, 9, 29}, 1, 5), // d = 4
        Parameters(
                // p, r, c, bits, h
                127, 2, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {686, 838}, {6, 17}, {9, 137}, 1, 5), // d = 4
        Parameters(
                // p, r, c, bits, h
                65537, 1, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {923, 475}, {2, 23}, {3, 461}, 1, 5), // d = 3
        Parameters(
                // p, r, c, bits, h
                31, 1, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {332, 1331}, {4, 110}, {5, 331}, 1, 5), // d = 3
        // idx = 4, p = 1 mod 4
        GTestFatboot::genPo2Param(10, 0, 1, 2, force_chen_han),
        // idx = 5, p = 3 mod 4
        GTestFatboot::genPo2Param(10, 2, -1, 2, force_chen_han),
        // thib bits checklist:
        //  1. p % 4 = 1: ok
        //  2. p % 4 = 3 & not forceRadix2: ok
        //  3. p % 4 = 3 & forceRadix2: ok
    };
    Parameters eg_params[] = {
        Parameters( // M=32551
                // p, r, c, bits, h
                127, 2, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {7571, 28768}, {42, -54}, {43, 757}, 1, 5), // d = 3
        Parameters( // M=45551
                // p, r, c, bits, h
                17, 4, 3, 800, 0, 
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {19394, 37270}, {100, 10}, {101, 11, 41}, 1, 5), // d = 3
    };
    Parameters po2_real_params[] = {
        // 2^16, p = 1 mod 4, d = 1
        /**
         * extration time ~= 5, bits = 266
         * partition:     firstMap (time & bits), secondMap (time & bits)
         * ### thin
         * [0, 8, 16]         7.2+6.6 / 60        17.1+20.0  / 98
         * [0, 6, 11, 16]     1.7+0.7+1.6 / 79    4.3+2.4+4.8 / 134
         * [0, 6, 12, 16]     3.8=0.9+1.2+1.6 / 79    10.8=4.3+3.8+2.6 / 134  <<--
         * [0, 5, 8, 12, 16]  3.36 / 99           9.11  / 170
         * ### thick
         * [0, 6, 12, 16]     10.8 / 126          10.7 / 126
        */
        Parameters( // d = 1, m = 65536
                // t = 1, log(qks*R) ~= 71.6, encapHwt = 26, seclvl = 134.4
                // main seclvl = 81.13, total bits <= 1350
                // p, r, c, bits, h
                65537, 1, 3, 980, 0, //  bits = 1332
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {-1, 5}, {2, 16384}, {65536}, 1, 5),
        // 2^16, p = 3 mod 4, d = 8
        /**
         * extraction time ~= 5.5s, bits = 232
         * ### thin (non-forceRadix2)
         * [0, 8, 13]  6.8=2.7+4.1 / 52   24.4=13.5+10 / 88
         * [0, 7, 13]  5.4=3.2+2.1 / 53   19.7=7.1+11.6 / 89
         * [0, 5, 9, 13] 4.2=2.5+1.0+0.7 / 71   16.4=2.4+3.6+9.5 / 121
         * [0, 6, 10, 13] 3.5=1.2+1.0+1.2 / 70  13.1=4.0+3.3+4.9 / 120 <<--
         * 
         * ### thin (forceRadix2)
         * [0, 6, 10, 13] 2.9=0.6+1.1+1.2 / 68  10.3=3.8+3.8+2.2 / 116
         * [0, 5, 9, 13]  2.8=1.0+1.1+0.2 / 69  10.1=2.3+3.8+3.5 / 118 <<--
         * [0, 6, 9, 13]  2.9 / 10.6
         * 
         * ### thick (non-forceRadix2)
         * [0, 6, 10, 13] 9.4+3.3+4.8 / 148 
         * [0, 5, 10, 13] 6.4+4.9+4.8 / 148 <<--
         * 
         * ### thick (forceRadix2)
         * [0, 5, 10, 13] 23.3=6.2+9.8+5.3 / 157
         * [0, 5, 9, 13]  21.2=6.2+5.6+7.4 / 157 <<--
        */
        Parameters( // d = 8, m = 65536
                // t = 1, log(qks*R) ~= 63.60, encapHwt = 24, seclvl = 129.8
                // main seclvl = 81.13, total bits <= 1350
                // p, r, c, bits, h
                8191, 1, 3, 980, 0, //  bits = 1332
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {5}, {4096}, {65536}, 1, 5),
        // 2^16, p = 3 mod 4, d = 2
        /**
         * ### thin (non-forceRadix2)
         * [0, 5, 9, 15]: 8.9=6.8+1.7+0.5  24.4=1.6+4.6+17.9
         * [0, 6, 10, 15]: 8.4=6.5+1.0+0.8 23.0=2.9+2.9+17.0
         * [0, 7, 11, 15]: 5.5=3.2+0.7+1.5 16.7=5.5+2.5+8.5
         * [0, 7, 12, 15]: 4.2=1.6+1.2+1.5 13.5=5.3+3.8+4.2 <<--
         * ### thin (forceRadix2)
         * [0, 6, 10, 15]: 3.4=1.6+1.0+0.9 12.1=2.9+3.3+5.8 <<--
         * [0, 7, 11, 15]: 3.4=0.8+1.0+1.6 12.1=5.5+3.5+3.1
         * 
         * ### thick (non-forceRadix2)
         * [0, 7, 12, 15]: 13.8=5.4+3.8+4.3 <<--
         * ### thick (forceRadix2)
         * [0, 6, 10, 15]:
        */
        Parameters( // d = 2, m = 65536
                // t = 1, log(qks*R) ~= 73.75, encapHwt = 26, seclvl = 133.81
                // main seclvl = 81.13, total bits <= 1350
                // p, r, c, bits, h
                131071, 1, 3, 980, 0, //  bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {5}, {16384}, {65536}, 1, 5),
    };
    // log2p < 60
    // log2p = 22, 29, 42, 54
    // NOTE: KSS24---M=2^18, logq <= 3440, total logq = 3300 bits
    Parameters KSS24_params[] = {
        Parameters( // d = 1, m = 2^18, 22 bits of p
                // t = 1, log(qks*R) ~= 91.41, encapHwt = 20, seclvl = 135
                // main seclvl > 129, total bits = 3355
                // p, r, c, bits, h
                5767169L, 1, 3, 2500, 240, //  bits
                // nthreads,seed,useCache,c_m,
                    100, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {-1, 5}, {2, 65536}, {262144}, 1, 5),
        Parameters( // d = 1, m = 2^18, 29 bits of p
                // t = 1, log(qks*R) ~= 108.79, encapHwt = 20, seclvl = 132
                // main seclvl > 129 , total bits = 3355
                // p, r, c, bits, h
                529268737L, 1, 3, 2500, 240, //  bits
                // nthreads,seed,useCache,c_m,
                    100, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {-1, 5}, {2, 65536}, {262144}, 1, 5),
        Parameters( // d = 1, m = 2^18, 42 bits of p
                // t = 1, log(qks*R) ~= 93.85, encapHwt = 20, seclvl = 135
                // main seclvl > 129 , total bits <= 3312
                // p, r, c, bits, h
                4398052540417L, 1, 3, 2500, 240, //  bits
                // nthreads,seed,useCache,c_m,
                    100, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {-1, 5}, {2, 65536}, {262144}, 1, 5),
        Parameters( // d = 1, m = 2^18, 54 bits of p
                // t = 1, log(qks*R) ~= 109.85, encapHwt = 20, seclvl = 132
                // main seclvl > 129 , total bits <= 3298
                // p, r, c, bits, h
                18014398512365569L, 1, 3, 2500, 240, //  bits
                // nthreads,seed,useCache,c_m,
                    100, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {-1, 5}, {2, 65536}, {262144}, 1, 5)
    };
    Parameters real_params[] = {
        // NOTE: parameters used in GIKV22
        // TODO: update bits
        // NOTE: we don't have the code of GIKV23... 
        //  if we use encapsulated sparse key, this is no way to compare with them
        // XXX: always set the mvec entries as primes, otherwise bound (5) may fail
        // Parameters( // d = 40, C = 149, good dim, m = 45551
        //         // log(qks*R) ~= 72, encapHwt = 22(12), seclvl = 129.1(82.3)
        //         // main seclvl = 78, total bits = (1593, 84.01793423200083)
        //         // p, r, c, bits, h
        //         17, 4, 3, 1100, 0, 
        //         // nthreads,seed,useCache,c_m,
        //             1, 0, 1, 100, 
        //         // force_bsgs,force_hoist,chen_han,debug,scale
        //             0, 0, force_chen_han, 0, 0, 
        //         // global_gens,global_ords,global_mvec,outer_rep,inner_rep
        //         {19394, 37270}, {100, 10}, {101, 11, 41}, 1, 5),
        // Parameters( // d = 14, C = 161, bad dim, m = 32551
        //         // log(qks*R) ~= 66, encapHwt = 24(10), seclvl = 134.7(70.9)
        //         // main seclvl = 66, total bits = (1593, 69.76804235959568)
        //         // p, r, c, bits, h
        //         127, 2, 3, 1100, 0, 
        //         // nthreads,seed,useCache,c_m,
        //             1, 0, 1, 100, 
        //         // force_bsgs,force_hoist,chen_han,debug,scale
        //             0, 0, force_chen_han, 0, 0, 
        //         // global_gens,global_ords,global_mvec,outer_rep,inner_rep
        //         {7571, 28768}, {42, -54}, {43, 757}, 1, 5),
        // NOTE: new params for p = 17 and 127
/*A*/   Parameters( // d = 24, C = 105, good dim, m = 38309
                // t = 2, log(qks*R) ~= 71, encapHwt = 14(24), seclvl = 89.8(133.1)
                // native, log(qks*R) ~= 65(64), encapHwt = 12(24), seclvl = 81.3(136.0)
                // main seclvl = 82.5, total bits <= 1500
                // p, r, c, bits, h
                17, 4, 3, 1070, 0, // 1462 bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {26421, 23781}, {28, 55}, {29, 1321}, 1, 5),
/*B*/   Parameters( // d = 45, C = 113, good dim, m = 56647
                // t = 1, log(qks*R) ~= 67, encapHwt = 12(22), seclvl = 85.7(135.4)
                // native, log(qks*R) ~= 61(70), encapHwt = 12(22), seclvl = 87.3(134.1)
                // TODO: larger t's
                // main seclvl = 82.6, total bits <= 2250
                // p, r, c, bits, h
                127, 2, 3, 1630, 0, // 2253 bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {12249, 49026}, {36, 34}, {37, 1531}, 1, 5),
        // NOTE: parameters with larger p
/*C*/   Parameters( // d = 28, C = 159, bad dim, m = 55427
                // t = 1, log(qks*R) ~= 71, encapHwt = 12(22), seclvl = 85.5(133.4)
                // native, log(qks*R) ~= 65(65), encapHwt = 12(22), seclvl = 87.4(135.5)
                // main seclvl = 82.8, total bits = 2200
                // p, r, c, bits, h
                257, 2, 3, 1610, 0,  // 2176 bits
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {52850, 2581}, {42, -46}, {43, 1289}, 1, 5),
/*D*/   Parameters( // d = 14, C = 129, good dim, m = 45193
                // t = 1, log(qks*R) ~= 73, encapHwt = 12(24), seclvl = 81.9(136.7)
                // t = -1, log(qks*R) ~= 63, encapHwt = 12(22), seclvl = 83.8(131.8)
                // native, log(qks*R) ~= 67(67), encapHwt = 12(22), seclvl = 83.3(130.2)
                // main seclvl = 82.3, total bits = 1800
                // p, r, c, bits, h
                8191, 1, 3, 1320, 0, // 1803 bits for t = 1, 1787 bits for t = -1
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {26276, 32595}, {42, 75}, {43, 1051}, 1, 5),
/*E*/   Parameters( // d = 18, C = 141, good dim, m = 50731
                // t = 1, log(qks*R) ~= 82, encapHwt = 12(24), seclvl = 82.3(136.8)
                // t = -1, log(qks*R) ~= 68, encapHwt = 12(22), seclvl = 84.8(132.8)
                // native, log(qks*R) ~= 76(75), encapHwt = 12(22), seclvl = 83.0(130.6)
                // main seclvl = 82.3, total bits = 2050
                // p, r, c, bits, h
                65537, 1, 3, 1500, 0, // 2036 bits for t = 1, 2018 for t = -1
                // nthreads,seed,useCache,c_m,
                    1, 0, 1, 100, 
                // force_bsgs,force_hoist,chen_han,debug,scale
                    0, 0, force_chen_han, 0, 0, 
                // global_gens,global_ords,global_mvec,outer_rep,inner_rep
                {48117, 5239}, {96, 29}, {97, 523}, 1, 5),

        };
    // clang-format on
    bool noPrint = false;
    bool dry = false; // don't set to true...
    long encapSkHwt = h_arg;
    long btsC = 3;
    long t = t_arg; // -1 for non-power-of-p aux, 0 for auto-decided power-of-p aux, >0 for aux = p^t
    double btsScale = 8;
    bool newBtsFlag = newBts_flag;
    bool newKSFlag = newKS_flag;
    bool s2cFirstFlag = s2cFirst_flag;
    bool isThick = thick_flag;
    // for ordinary key with 128 bit security, use h = 492 (no, achieving 128 is too expensive)
    // use h = 120 for the ordinary key
    std::cout << "chosen idx = " << i_arg << "\n";
    std::cout << "encapHwt = " << encapSkHwt << "\n";
    std::cout << "t for new bts is " << t << "\n";
    std::cout << "new bts ? " << newBtsFlag << "\n";
    std::cout << "new ks ? " << newKSFlag << "\n";
    std::cout << "thick bts ?" << isThick << "\n";
    std::cout << "s2cFirst ? " << s2cFirstFlag << "\n";
    std::cout << "repeat times = " << repeat_arg << "\n";

	Parameters *dst_params = po2_real_params;
    // dst_params = KSS24_params;
    dst_params[i_arg].inner_rep = repeat_arg;
    dst_params[i_arg].useCache = cache_arg;
    dummy_arg = !thick_flag && baseline_arg;
    if (dummy_arg) {
        // dst_params[i_arg].useCache = false; // prevent OOM
        if (dst_params[i_arg].global_gens.length() != 2){
            long cur_ord = dst_params[i_arg].global_ords[0];
            long m = helib::computeProd(dst_params[i_arg].global_mvec);

            dst_params[i_arg].global_gens.SetLength(2);
            dst_params[i_arg].global_gens[0] = NTL::PowerMod(5, cur_ord/2, m);
            dst_params[i_arg].global_gens[1] = 5;
            dst_params[i_arg].global_ords.SetLength(2);
            dst_params[i_arg].global_ords[0] = 2;
            dst_params[i_arg].global_ords[1] = cur_ord / 2;
        }
    }
    std::vector<long> partition;
    long m = helib::computeProd(dst_params[i_arg].global_mvec);
    long p = dst_params[i_arg].p;
    if (helib::is2power(m)) {
        // partition = auto_po2_parition(m, p, forceRadix2_arg);

        // XXX: manual partition here
        if (dst_params == po2_real_params){
            if (i_arg == 0)
                partition = {0, 6, 12, 16};
            if (i_arg == 1) {
                if (!isThick && !forceRadix2_arg)
                    partition = {0, 6, 10, 13};
                else if (isThick && !forceRadix2_arg)
                    partition = {0, 5, 10, 13};
                else
                    partition = {0, 5, 9, 13};
            }
            if (i_arg == 2) {
                if (forceRadix2_arg)
                    partition = {0, 6, 10, 15};
                else
                    partition = {0, 7, 12, 15};
            }

            if (isThick && baseline_arg) {
                if (i_arg == 0)
                    partition = {0, 16};
                if (i_arg == 1)
                    partition = {0, 13};
                if (i_arg == 2)
                    partition = {0, 15};
            }
		}
        // XXX: KSS24, M=2^18
        // partition = {0, 7, 13, 18};

        std::cout << "### partition is ";
        for (long pt : partition)
            std::cout << pt << " ";
        std::cout << "\n";
    }

    GTestFatboot test(dst_params[i_arg], noPrint, dry, encapSkHwt, 
        btsC, t, btsScale, newBtsFlag, newKSFlag, s2cFirstFlag, isThick, 
        forceRadix2_arg, partition, dummy_arg);
    test.run();
    test.TearDown();
    return 0;
}

// LEGACY TEST DEFAULT PARAMETERS:
// long p=2;
// long r=1;
// long c=3;
// long bits=600;
// long t=64;
// long nthreads=1;
// long seed=0;
// long useCache=1;
// int c_m = 100;
