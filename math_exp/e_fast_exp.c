#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "e_lib.h"

// GCC
#define RESTRICT __restrict__

extern float expapprox(float val);
extern void expapprox4(float *RESTRICT dst, const float *src);

//
// -------------------------------------------------------------------------------------
// Fast approximate exp() function for Epiphany.
//
// Copyright 2015 Syoyo Fujita <syoyo@lighttransport.com>
// Licensed under Apache 2.0 License.
//
//
// <<Epiphany>>
//   -O3 -std=c99 -fsingle-precision-constant -mno-soft-cmpsf -mcmove
//   -mfp-mode=truncate
//
//                          Enable range check       Disable range check
//  -------------+--------+-------------------------+-------------------------------
//  expapprox()  | scalar | 74 cycles               | 54 cycles
//  expapprox4() | 4 SIMD | 32 cycles(128 in total) | 25 cycles(90 in total)
//
//
//  fmath_exp()  | scalar | 33 cycles
//  fmath_exp4() | 4 SIMD | 18 cycles(74 in total)
//
// 
//  fmath_exp() is faster and more accurate than expapprox(), but with the cost of
//  table buffer(512 byte ~ 4KB)
//
//  - Reference
//
//  -------------+--------+---------------------------------------------------------
//  expf()       | scalar | 141633 cycles(tooooo slow because of SW implementation?)
//
//  - Difference
//
//    * expapprox()
//      [-30.0, 30.0] Relative diff.
//
//      ave = 2.429884e-06, min = 0.000000e+00, max = 7.674494e-06
//
//    * fmath_exp() tablesize = 10(4KB)
//
//      [-30, 30]
//
//      table size    ave rel. diff    min rel. diff    max rel. diff
//      ------------+----------------+----------------+-----------------------
//      10 (4KB)    | 7.261362e-08   | 0.000000e+00   | 1.078697e-06
//      8  (1KB)    | 3.263584e-07   | 0.000000e+00   | 1.651098e-06
//      7  (512B)   | 1.233390e-06   | 0.000000e+00   | 3.771282e-06
//
//  - Note
//    * -mfp-mode=truncate may loose some precison, but emits more optimal
//    assembly.
//    * -mno-soft-compsf required to not emit software compare function, which
//    breaks
//      some IEEE 754 compliance, but should be totally OK for approximate math
//      functions.
//
//

// If you are sure that the input value for exp() is within [-88.0, 88.0],
// you can disable range check, which results in faster evaluation.
#define FMATH_EXP_DISABLE_RANGE_CHECK (0)

// base on fmath::exp https://github.com/herumi/fmath/blob/master/fmath.hpp

#define FMATH_EXP_TABLE_SIZE	(7)

#if (FMATH_EXP_TABLE_SIZE == 10)
const unsigned int kFmathExpTable[1024] = {
  0x00000000, 0x00001630, 0x00002c64, 0x0000429c, 
  0x000058d8, 0x00006f17, 0x0000855b, 0x00009ba2, 
  0x0000b1ed, 0x0000c83c, 0x0000de8f, 0x0000f4e6, 
  0x00010b41, 0x0001219f, 0x00013802, 0x00014e68, 
  0x000164d2, 0x00017b40, 0x000191b2, 0x0001a828, 
  0x0001bea1, 0x0001d51f, 0x0001eba1, 0x00020226, 
  0x000218af, 0x00022f3c, 0x000245ce, 0x00025c63, 
  0x000272fc, 0x00028998, 0x0002a039, 0x0002b6de, 
  0x0002cd87, 0x0002e433, 0x0002fae4, 0x00031198, 
  0x00032850, 0x00033f0d, 0x000355cd, 0x00036c91, 
  0x00038359, 0x00039a25, 0x0003b0f5, 0x0003c7c9, 
  0x0003dea1, 0x0003f57d, 0x00040c5d, 0x00042341, 
  0x00043a29, 0x00045115, 0x00046804, 0x00047ef8, 
  0x000495f0, 0x0004aceb, 0x0004c3eb, 0x0004daef, 
  0x0004f1f6, 0x00050902, 0x00052012, 0x00053725, 
  0x00054e3d, 0x00056558, 0x00057c78, 0x0005939c, 
  0x0005aac3, 0x0005c1ef, 0x0005d91f, 0x0005f052, 
  0x0006078a, 0x00061ec6, 0x00063606, 0x00064d4a, 
  0x00066491, 0x00067bdd, 0x0006932d, 0x0006aa81, 
  0x0006c1d9, 0x0006d935, 0x0006f095, 0x000707f9, 
  0x00071f62, 0x000736ce, 0x00074e3e, 0x000765b3, 
  0x00077d2b, 0x000794a8, 0x0007ac28, 0x0007c3ad, 
  0x0007db35, 0x0007f2c2, 0x00080a53, 0x000821e8, 
  0x00083981, 0x0008511e, 0x000868c0, 0x00088065, 
  0x0008980f, 0x0008afbc, 0x0008c76e, 0x0008df23, 
  0x0008f6dd, 0x00090e9b, 0x0009265d, 0x00093e24, 
  0x000955ee, 0x00096dbc, 0x0009858f, 0x00099d66, 
  0x0009b541, 0x0009cd20, 0x0009e503, 0x0009fcea, 
  0x000a14d5, 0x000a2cc5, 0x000a44b9, 0x000a5cb1, 
  0x000a74ad, 0x000a8cad, 0x000aa4b1, 0x000abcba, 
  0x000ad4c6, 0x000aecd7, 0x000b04ec, 0x000b1d05, 
  0x000b3523, 0x000b4d44, 0x000b656a, 0x000b7d94, 
  0x000b95c2, 0x000badf4, 0x000bc62b, 0x000bde65, 
  0x000bf6a4, 0x000c0ee7, 0x000c272f, 0x000c3f7a, 
  0x000c57ca, 0x000c701e, 0x000c8876, 0x000ca0d2, 
  0x000cb933, 0x000cd198, 0x000cea01, 0x000d026e, 
  0x000d1adf, 0x000d3355, 0x000d4bcf, 0x000d644d, 
  0x000d7cd0, 0x000d9556, 0x000dade1, 0x000dc671, 
  0x000ddf04, 0x000df79c, 0x000e1038, 0x000e28d8, 
  0x000e417d, 0x000e5a25, 0x000e72d3, 0x000e8b84, 
  0x000ea43a, 0x000ebcf3, 0x000ed5b2, 0x000eee74, 
  0x000f073b, 0x000f2006, 0x000f38d5, 0x000f51a9, 
  0x000f6a81, 0x000f835d, 0x000f9c3e, 0x000fb523, 
  0x000fce0c, 0x000fe6fa, 0x000fffec, 0x001018e2, 
  0x001031dc, 0x00104adb, 0x001063de, 0x00107ce6, 
  0x001095f2, 0x0010af02, 0x0010c816, 0x0010e12f, 
  0x0010fa4d, 0x0011136e, 0x00112c94, 0x001145be, 
  0x00115eed, 0x00117820, 0x00119158, 0x0011aa93, 
  0x0011c3d3, 0x0011dd18, 0x0011f661, 0x00120fae, 
  0x00122900, 0x00124256, 0x00125bb0, 0x0012750f, 
  0x00128e72, 0x0012a7da, 0x0012c146, 0x0012dab7, 
  0x0012f42c, 0x00130da5, 0x00132723, 0x001340a5, 
  0x00135a2b, 0x001373b6, 0x00138d46, 0x0013a6d9, 
  0x0013c072, 0x0013da0e, 0x0013f3af, 0x00140d55, 
  0x001426ff, 0x001440ae, 0x00145a60, 0x00147418, 
  0x00148dd4, 0x0014a794, 0x0014c159, 0x0014db22, 
  0x0014f4f0, 0x00150ec2, 0x00152898, 0x00154274, 
  0x00155c53, 0x00157637, 0x00159020, 0x0015aa0d, 
  0x0015c3ff, 0x0015ddf5, 0x0015f7ef, 0x001611ee, 
  0x00162bf2, 0x001645fa, 0x00166006, 0x00167a18, 
  0x0016942d, 0x0016ae47, 0x0016c866, 0x0016e289, 
  0x0016fcb1, 0x001716dd, 0x0017310e, 0x00174b43, 
  0x0017657d, 0x00177fbc, 0x001799ff, 0x0017b446, 
  0x0017ce92, 0x0017e8e3, 0x00180338, 0x00181d92, 
  0x001837f0, 0x00185253, 0x00186cbb, 0x00188727, 
  0x0018a197, 0x0018bc0d, 0x0018d686, 0x0018f105, 
  0x00190b88, 0x0019260f, 0x0019409c, 0x00195b2c, 
  0x001975c2, 0x0019905c, 0x0019aafa, 0x0019c59e, 
  0x0019e046, 0x0019faf2, 0x001a15a3, 0x001a3059, 
  0x001a4b13, 0x001a65d2, 0x001a8096, 0x001a9b5e, 
  0x001ab62b, 0x001ad0fd, 0x001aebd3, 0x001b06ae, 
  0x001b218d, 0x001b3c71, 0x001b575a, 0x001b7248, 
  0x001b8d3a, 0x001ba831, 0x001bc32c, 0x001bde2c, 
  0x001bf931, 0x001c143b, 0x001c2f49, 0x001c4a5c, 
  0x001c6573, 0x001c8090, 0x001c9bb1, 0x001cb6d6, 
  0x001cd201, 0x001ced30, 0x001d0864, 0x001d239c, 
  0x001d3eda, 0x001d5a1c, 0x001d7562, 0x001d90ae, 
  0x001dabfe, 0x001dc753, 0x001de2ad, 0x001dfe0b, 
  0x001e196e, 0x001e34d6, 0x001e5043, 0x001e6bb4, 
  0x001e872a, 0x001ea2a5, 0x001ebe25, 0x001ed9a9, 
  0x001ef532, 0x001f10c0, 0x001f2c53, 0x001f47eb, 
  0x001f6387, 0x001f7f28, 0x001f9ace, 0x001fb679, 
  0x001fd228, 0x001feddc, 0x00200996, 0x00202553, 
  0x00204116, 0x00205cde, 0x002078aa, 0x0020947b, 
  0x0020b051, 0x0020cc2c, 0x0020e80b, 0x002103f0, 
  0x00211fd9, 0x00213bc7, 0x002157ba, 0x002173b2, 
  0x00218faf, 0x0021abb0, 0x0021c7b7, 0x0021e3c2, 
  0x0021ffd2, 0x00221be7, 0x00223801, 0x0022541f, 
  0x00227043, 0x00228c6b, 0x0022a899, 0x0022c4cb, 
  0x0022e102, 0x0022fd3e, 0x0023197f, 0x002335c5, 
  0x0023520f, 0x00236e5f, 0x00238ab3, 0x0023a70d, 
  0x0023c36b, 0x0023dfce, 0x0023fc37, 0x002418a4, 
  0x00243516, 0x0024518d, 0x00246e08, 0x00248a89, 
  0x0024a70f, 0x0024c39a, 0x0024e029, 0x0024fcbe, 
  0x00251958, 0x002535f6, 0x00255299, 0x00256f42, 
  0x00258bef, 0x0025a8a2, 0x0025c559, 0x0025e215, 
  0x0025fed7, 0x00261b9d, 0x00263868, 0x00265538, 
  0x0026720e, 0x00268ee8, 0x0026abc7, 0x0026c8ac, 
  0x0026e595, 0x00270283, 0x00271f76, 0x00273c6f, 
  0x0027596c, 0x0027766e, 0x00279376, 0x0027b082, 
  0x0027cd94, 0x0027eaaa, 0x002807c6, 0x002824e6, 
  0x0028420c, 0x00285f37, 0x00287c66, 0x0028999b, 
  0x0028b6d5, 0x0028d414, 0x0028f158, 0x00290ea1, 
  0x00292bef, 0x00294942, 0x0029669b, 0x002983f8, 
  0x0029a15b, 0x0029bec2, 0x0029dc2f, 0x0029f9a1, 
  0x002a1718, 0x002a3494, 0x002a5215, 0x002a6f9b, 
  0x002a8d26, 0x002aaab7, 0x002ac84c, 0x002ae5e7, 
  0x002b0387, 0x002b212c, 0x002b3ed6, 0x002b5c85, 
  0x002b7a3a, 0x002b97f3, 0x002bb5b2, 0x002bd376, 
  0x002bf13f, 0x002c0f0d, 0x002c2ce0, 0x002c4ab9, 
  0x002c6897, 0x002c867a, 0x002ca462, 0x002cc24f, 
  0x002ce041, 0x002cfe39, 0x002d1c36, 0x002d3a38, 
  0x002d583f, 0x002d764b, 0x002d945d, 0x002db274, 
  0x002dd090, 0x002deeb1, 0x002e0cd8, 0x002e2b03, 
  0x002e4934, 0x002e676b, 0x002e85a6, 0x002ea3e7, 
  0x002ec22d, 0x002ee078, 0x002efec8, 0x002f1d1e, 
  0x002f3b79, 0x002f59d9, 0x002f783e, 0x002f96a9, 
  0x002fb519, 0x002fd38e, 0x002ff209, 0x00301089, 
  0x00302f0e, 0x00304d98, 0x00306c28, 0x00308abd, 
  0x0030a957, 0x0030c7f7, 0x0030e69c, 0x00310546, 
  0x003123f6, 0x003142aa, 0x00316165, 0x00318024, 
  0x00319ee9, 0x0031bdb3, 0x0031dc83, 0x0031fb57, 
  0x00321a32, 0x00323911, 0x003257f6, 0x003276e0, 
  0x003295d0, 0x0032b4c5, 0x0032d3bf, 0x0032f2bf, 
  0x003311c4, 0x003330cf, 0x00334fde, 0x00336ef4, 
  0x00338e0e, 0x0033ad2e, 0x0033cc54, 0x0033eb7e, 
  0x00340aaf, 0x003429e4, 0x0034491f, 0x00346860, 
  0x003487a6, 0x0034a6f1, 0x0034c642, 0x0034e598, 
  0x003504f3, 0x00352454, 0x003543bb, 0x00356327, 
  0x00358298, 0x0035a20f, 0x0035c18b, 0x0035e10d, 
  0x00360094, 0x00362020, 0x00363fb2, 0x00365f4a, 
  0x00367ee7, 0x00369e89, 0x0036be31, 0x0036dddf, 
  0x0036fd92, 0x00371d4a, 0x00373d08, 0x00375ccc, 
  0x00377c95, 0x00379c63, 0x0037bc37, 0x0037dc11, 
  0x0037fbf0, 0x00381bd4, 0x00383bbe, 0x00385bae, 
  0x00387ba3, 0x00389b9e, 0x0038bb9e, 0x0038dba4, 
  0x0038fbaf, 0x00391bc0, 0x00393bd7, 0x00395bf3, 
  0x00397c14, 0x00399c3b, 0x0039bc68, 0x0039dc9a, 
  0x0039fcd2, 0x003a1d10, 0x003a3d53, 0x003a5d9b, 
  0x003a7dea, 0x003a9e3e, 0x003abe97, 0x003adef6, 
  0x003aff5b, 0x003b1fc5, 0x003b4035, 0x003b60aa, 
  0x003b8126, 0x003ba1a6, 0x003bc22d, 0x003be2b9, 
  0x003c034a, 0x003c23e2, 0x003c447f, 0x003c6521, 
  0x003c85ca, 0x003ca678, 0x003cc72b, 0x003ce7e5, 
  0x003d08a4, 0x003d2968, 0x003d4a33, 0x003d6b03, 
  0x003d8bd8, 0x003dacb4, 0x003dcd95, 0x003dee7c, 
  0x003e0f68, 0x003e305a, 0x003e5152, 0x003e7250, 
  0x003e9353, 0x003eb45c, 0x003ed56b, 0x003ef67f, 
  0x003f179a, 0x003f38ba, 0x003f59df, 0x003f7b0b, 
  0x003f9c3c, 0x003fbd73, 0x003fdeb0, 0x003ffff2, 
  0x0040213b, 0x00404289, 0x004063dc, 0x00408536, 
  0x0040a695, 0x0040c7fb, 0x0040e966, 0x00410ad6, 
  0x00412c4d, 0x00414dc9, 0x00416f4b, 0x004190d3, 
  0x0041b261, 0x0041d3f5, 0x0041f58e, 0x0042172d, 
  0x004238d2, 0x00425a7d, 0x00427c2e, 0x00429de4, 
  0x0042bfa1, 0x0042e163, 0x0043032b, 0x004324f9, 
  0x004346cd, 0x004368a7, 0x00438a86, 0x0043ac6b, 
  0x0043ce57, 0x0043f048, 0x0044123f, 0x0044343c, 
  0x0044563f, 0x00447848, 0x00449a56, 0x0044bc6b, 
  0x0044de85, 0x004500a5, 0x004522cc, 0x004544f8, 
  0x0045672a, 0x00458962, 0x0045aba0, 0x0045cde4, 
  0x0045f02e, 0x0046127e, 0x004634d3, 0x0046572f, 
  0x00467991, 0x00469bf8, 0x0046be66, 0x0046e0d9, 
  0x00470353, 0x004725d2, 0x00474858, 0x00476ae3, 
  0x00478d75, 0x0047b00c, 0x0047d2aa, 0x0047f54d, 
  0x004817f7, 0x00483aa6, 0x00485d5b, 0x00488017, 
  0x0048a2d8, 0x0048c5a0, 0x0048e86d, 0x00490b41, 
  0x00492e1b, 0x004950fa, 0x004973e0, 0x004996cc, 
  0x0049b9be, 0x0049dcb5, 0x0049ffb3, 0x004a22b7, 
  0x004a45c1, 0x004a68d1, 0x004a8be8, 0x004aaf04, 
  0x004ad226, 0x004af54f, 0x004b187d, 0x004b3bb2, 
  0x004b5eed, 0x004b822e, 0x004ba575, 0x004bc8c2, 
  0x004bec15, 0x004c0f6e, 0x004c32ce, 0x004c5633, 
  0x004c799f, 0x004c9d11, 0x004cc089, 0x004ce407, 
  0x004d078c, 0x004d2b16, 0x004d4ea7, 0x004d723d, 
  0x004d95da, 0x004db97e, 0x004ddd27, 0x004e00d6, 
  0x004e248c, 0x004e4848, 0x004e6c0a, 0x004e8fd2, 
  0x004eb3a1, 0x004ed775, 0x004efb50, 0x004f1f31, 
  0x004f4319, 0x004f6706, 0x004f8afa, 0x004faef4, 
  0x004fd2f4, 0x004ff6fb, 0x00501b08, 0x00503f1b, 
  0x00506334, 0x00508753, 0x0050ab79, 0x0050cfa5, 
  0x0050f3d7, 0x00511810, 0x00513c4f, 0x00516094, 
  0x005184df, 0x0051a931, 0x0051cd89, 0x0051f1e7, 
  0x0052164c, 0x00523ab7, 0x00525f28, 0x005283a0, 
  0x0052a81e, 0x0052cca2, 0x0052f12c, 0x005315bd, 
  0x00533a54, 0x00535ef2, 0x00538396, 0x0053a840, 
  0x0053ccf1, 0x0053f1a8, 0x00541665, 0x00543b29, 
  0x00545ff3, 0x005484c3, 0x0054a99a, 0x0054ce77, 
  0x0054f35b, 0x00551845, 0x00553d35, 0x0055622c, 
  0x00558729, 0x0055ac2d, 0x0055d137, 0x0055f647, 
  0x00561b5e, 0x0056407b, 0x0056659f, 0x00568ac9, 
  0x0056affa, 0x0056d531, 0x0056fa6e, 0x00571fb2, 
  0x005744fd, 0x00576a4e, 0x00578fa5, 0x0057b503, 
  0x0057da67, 0x0057ffd2, 0x00582543, 0x00584abb, 
  0x00587039, 0x005895be, 0x0058bb49, 0x0058e0db, 
  0x00590673, 0x00592c12, 0x005951b8, 0x00597763, 
  0x00599d16, 0x0059c2cf, 0x0059e88e, 0x005a0e54, 
  0x005a3421, 0x005a59f4, 0x005a7fcd, 0x005aa5ae, 
  0x005acb94, 0x005af182, 0x005b1776, 0x005b3d70, 
  0x005b6371, 0x005b8979, 0x005baf87, 0x005bd59c, 
  0x005bfbb8, 0x005c21da, 0x005c4802, 0x005c6e32, 
  0x005c9468, 0x005cbaa4, 0x005ce0e7, 0x005d0731, 
  0x005d2d82, 0x005d53d9, 0x005d7a36, 0x005da09b, 
  0x005dc706, 0x005ded77, 0x005e13f0, 0x005e3a6f, 
  0x005e60f5, 0x005e8781, 0x005eae14, 0x005ed4ae, 
  0x005efb4e, 0x005f21f5, 0x005f48a3, 0x005f6f58, 
  0x005f9613, 0x005fbcd5, 0x005fe39e, 0x00600a6d, 
  0x00603143, 0x00605820, 0x00607f03, 0x0060a5ee, 
  0x0060ccdf, 0x0060f3d7, 0x00611ad5, 0x006141db, 
  0x006168e7, 0x00618ffa, 0x0061b713, 0x0061de34, 
  0x0062055b, 0x00622c89, 0x006253be, 0x00627af9, 
  0x0062a23c, 0x0062c985, 0x0062f0d5, 0x0063182c, 
  0x00633f89, 0x006366ee, 0x00638e59, 0x0063b5cb, 
  0x0063dd44, 0x006404c4, 0x00642c4b, 0x006453d8, 
  0x00647b6d, 0x0064a308, 0x0064caaa, 0x0064f253, 
  0x00651a03, 0x006541b9, 0x00656977, 0x0065913c, 
  0x0065b907, 0x0065e0d9, 0x006608b2, 0x00663092, 
  0x00665879, 0x00668067, 0x0066a85c, 0x0066d058, 
  0x0066f85b, 0x00672064, 0x00674875, 0x0067708c, 
  0x006798ab, 0x0067c0d0, 0x0067e8fd, 0x00681130, 
  0x0068396a, 0x006861ac, 0x006889f4, 0x0068b243, 
  0x0068da99, 0x006902f7, 0x00692b5b, 0x006953c6, 
  0x00697c38, 0x0069a4b1, 0x0069cd32, 0x0069f5b9, 
  0x006a1e47, 0x006a46dd, 0x006a6f79, 0x006a981c, 
  0x006ac0c7, 0x006ae978, 0x006b1231, 0x006b3af1, 
  0x006b63b7, 0x006b8c85, 0x006bb55a, 0x006bde36, 
  0x006c0719, 0x006c3003, 0x006c58f4, 0x006c81ec, 
  0x006caaec, 0x006cd3f2, 0x006cfd00, 0x006d2614, 
  0x006d4f30, 0x006d7853, 0x006da17d, 0x006dcaae, 
  0x006df3e7, 0x006e1d26, 0x006e466d, 0x006e6fbb, 
  0x006e9910, 0x006ec26c, 0x006eebcf, 0x006f1539, 
  0x006f3eab, 0x006f6824, 0x006f91a4, 0x006fbb2b, 
  0x006fe4ba, 0x00700e4f, 0x007037ec, 0x00706190, 
  0x00708b3b, 0x0070b4ee, 0x0070dea8, 0x00710868, 
  0x00713231, 0x00715c00, 0x007185d7, 0x0071afb5, 
  0x0071d99a, 0x00720386, 0x00722d7a, 0x00725775, 
  0x00728177, 0x0072ab81, 0x0072d592, 0x0072ffaa, 
  0x007329c9, 0x007353f0, 0x00737e1e, 0x0073a853, 
  0x0073d290, 0x0073fcd4, 0x0074271f, 0x00745172, 
  0x00747bcc, 0x0074a62d, 0x0074d096, 0x0074fb06, 
  0x0075257d, 0x00754ffc, 0x00757a82, 0x0075a50f, 
  0x0075cfa4, 0x0075fa40, 0x007624e4, 0x00764f8f, 
  0x00767a41, 0x0076a4fb, 0x0076cfbc, 0x0076fa85, 
  0x00772555, 0x0077502d, 0x00777b0b, 0x0077a5f2, 
  0x0077d0df, 0x0077fbd5, 0x007826d1, 0x007851d5, 
  0x00787ce1, 0x0078a7f4, 0x0078d30e, 0x0078fe30, 
  0x0079295a, 0x0079548b, 0x00797fc3, 0x0079ab03, 
  0x0079d64a, 0x007a0199, 0x007a2cf0, 0x007a584d, 
  0x007a83b3, 0x007aaf20, 0x007ada94, 0x007b0610, 
  0x007b3194, 0x007b5d1f, 0x007b88b2, 0x007bb44c, 
  0x007bdfed, 0x007c0b97, 0x007c3748, 0x007c6300, 
  0x007c8ec0, 0x007cba88, 0x007ce657, 0x007d122e, 
  0x007d3e0c, 0x007d69f2, 0x007d95e0, 0x007dc1d5, 
  0x007dedd2, 0x007e19d6, 0x007e45e2, 0x007e71f6, 
  0x007e9e11, 0x007eca34, 0x007ef65f, 0x007f2291, 
  0x007f4ecb, 0x007f7b0d, 0x007fa756, 0x007fd3a7
};
#elif (FMATH_EXP_TABLE_SIZE == 8)
const unsigned int kFmathExpTable[256] = {
  0x00000000, 0x000058d8, 0x0000b1ed, 0x00010b41, 
  0x000164d2, 0x0001bea1, 0x000218af, 0x000272fc, 
  0x0002cd87, 0x00032850, 0x00038359, 0x0003dea1, 
  0x00043a29, 0x000495f0, 0x0004f1f6, 0x00054e3d, 
  0x0005aac3, 0x0006078a, 0x00066491, 0x0006c1d9, 
  0x00071f62, 0x00077d2b, 0x0007db35, 0x00083981, 
  0x0008980f, 0x0008f6dd, 0x000955ee, 0x0009b541, 
  0x000a14d5, 0x000a74ad, 0x000ad4c6, 0x000b3523, 
  0x000b95c2, 0x000bf6a4, 0x000c57ca, 0x000cb933, 
  0x000d1adf, 0x000d7cd0, 0x000ddf04, 0x000e417d, 
  0x000ea43a, 0x000f073b, 0x000f6a81, 0x000fce0c, 
  0x001031dc, 0x001095f2, 0x0010fa4d, 0x00115eed, 
  0x0011c3d3, 0x00122900, 0x00128e72, 0x0012f42c, 
  0x00135a2b, 0x0013c072, 0x001426ff, 0x00148dd4, 
  0x0014f4f0, 0x00155c53, 0x0015c3ff, 0x00162bf2, 
  0x0016942d, 0x0016fcb1, 0x0017657d, 0x0017ce92, 
  0x001837f0, 0x0018a197, 0x00190b88, 0x001975c2, 
  0x0019e046, 0x001a4b13, 0x001ab62b, 0x001b218d, 
  0x001b8d3a, 0x001bf931, 0x001c6573, 0x001cd201, 
  0x001d3eda, 0x001dabfe, 0x001e196e, 0x001e872a, 
  0x001ef532, 0x001f6387, 0x001fd228, 0x00204116, 
  0x0020b051, 0x00211fd9, 0x00218faf, 0x0021ffd2, 
  0x00227043, 0x0022e102, 0x0023520f, 0x0023c36b, 
  0x00243516, 0x0024a70f, 0x00251958, 0x00258bef, 
  0x0025fed7, 0x0026720e, 0x0026e595, 0x0027596c, 
  0x0027cd94, 0x0028420c, 0x0028b6d5, 0x00292bef, 
  0x0029a15b, 0x002a1718, 0x002a8d26, 0x002b0387, 
  0x002b7a3a, 0x002bf13f, 0x002c6897, 0x002ce041, 
  0x002d583f, 0x002dd090, 0x002e4934, 0x002ec22d, 
  0x002f3b79, 0x002fb519, 0x00302f0e, 0x0030a957, 
  0x003123f6, 0x00319ee9, 0x00321a32, 0x003295d0, 
  0x003311c4, 0x00338e0e, 0x00340aaf, 0x003487a6, 
  0x003504f3, 0x00358298, 0x00360094, 0x00367ee7, 
  0x0036fd92, 0x00377c95, 0x0037fbf0, 0x00387ba3, 
  0x0038fbaf, 0x00397c14, 0x0039fcd2, 0x003a7dea, 
  0x003aff5b, 0x003b8126, 0x003c034a, 0x003c85ca, 
  0x003d08a4, 0x003d8bd8, 0x003e0f68, 0x003e9353, 
  0x003f179a, 0x003f9c3c, 0x0040213b, 0x0040a695, 
  0x00412c4d, 0x0041b261, 0x004238d2, 0x0042bfa1, 
  0x004346cd, 0x0043ce57, 0x0044563f, 0x0044de85, 
  0x0045672a, 0x0045f02e, 0x00467991, 0x00470353, 
  0x00478d75, 0x004817f7, 0x0048a2d8, 0x00492e1b, 
  0x0049b9be, 0x004a45c1, 0x004ad226, 0x004b5eed, 
  0x004bec15, 0x004c799f, 0x004d078c, 0x004d95da, 
  0x004e248c, 0x004eb3a1, 0x004f4319, 0x004fd2f4, 
  0x00506334, 0x0050f3d7, 0x005184df, 0x0052164c, 
  0x0052a81e, 0x00533a54, 0x0053ccf1, 0x00545ff3, 
  0x0054f35b, 0x00558729, 0x00561b5e, 0x0056affa, 
  0x005744fd, 0x0057da67, 0x00587039, 0x00590673, 
  0x00599d16, 0x005a3421, 0x005acb94, 0x005b6371, 
  0x005bfbb8, 0x005c9468, 0x005d2d82, 0x005dc706, 
  0x005e60f5, 0x005efb4e, 0x005f9613, 0x00603143, 
  0x0060ccdf, 0x006168e7, 0x0062055b, 0x0062a23c, 
  0x00633f89, 0x0063dd44, 0x00647b6d, 0x00651a03, 
  0x0065b907, 0x00665879, 0x0066f85b, 0x006798ab, 
  0x0068396a, 0x0068da99, 0x00697c38, 0x006a1e47, 
  0x006ac0c7, 0x006b63b7, 0x006c0719, 0x006caaec, 
  0x006d4f30, 0x006df3e7, 0x006e9910, 0x006f3eab, 
  0x006fe4ba, 0x00708b3b, 0x00713231, 0x0071d99a, 
  0x00728177, 0x007329c9, 0x0073d290, 0x00747bcc, 
  0x0075257d, 0x0075cfa4, 0x00767a41, 0x00772555, 
  0x0077d0df, 0x00787ce1, 0x0079295a, 0x0079d64a, 
  0x007a83b3, 0x007b3194, 0x007bdfed, 0x007c8ec0, 
  0x007d3e0c, 0x007dedd2, 0x007e9e11, 0x007f4ecb
};
#elif (FMATH_EXP_TABLE_SIZE == 7)
const unsigned int kFmathExpTable[128] = {
  0x00000000, 0x0000b1ed, 0x000164d2, 0x000218af, 
  0x0002cd87, 0x00038359, 0x00043a29, 0x0004f1f6, 
  0x0005aac3, 0x00066491, 0x00071f62, 0x0007db35, 
  0x0008980f, 0x000955ee, 0x000a14d5, 0x000ad4c6, 
  0x000b95c2, 0x000c57ca, 0x000d1adf, 0x000ddf04, 
  0x000ea43a, 0x000f6a81, 0x001031dc, 0x0010fa4d, 
  0x0011c3d3, 0x00128e72, 0x00135a2b, 0x001426ff, 
  0x0014f4f0, 0x0015c3ff, 0x0016942d, 0x0017657d, 
  0x001837f0, 0x00190b88, 0x0019e046, 0x001ab62b, 
  0x001b8d3a, 0x001c6573, 0x001d3eda, 0x001e196e, 
  0x001ef532, 0x001fd228, 0x0020b051, 0x00218faf, 
  0x00227043, 0x0023520f, 0x00243516, 0x00251958, 
  0x0025fed7, 0x0026e595, 0x0027cd94, 0x0028b6d5, 
  0x0029a15b, 0x002a8d26, 0x002b7a3a, 0x002c6897, 
  0x002d583f, 0x002e4934, 0x002f3b79, 0x00302f0e, 
  0x003123f6, 0x00321a32, 0x003311c4, 0x00340aaf, 
  0x003504f3, 0x00360094, 0x0036fd92, 0x0037fbf0, 
  0x0038fbaf, 0x0039fcd2, 0x003aff5b, 0x003c034a, 
  0x003d08a4, 0x003e0f68, 0x003f179a, 0x0040213b, 
  0x00412c4d, 0x004238d2, 0x004346cd, 0x0044563f, 
  0x0045672a, 0x00467991, 0x00478d75, 0x0048a2d8, 
  0x0049b9be, 0x004ad226, 0x004bec15, 0x004d078c, 
  0x004e248c, 0x004f4319, 0x00506334, 0x005184df, 
  0x0052a81e, 0x0053ccf1, 0x0054f35b, 0x00561b5e, 
  0x005744fd, 0x00587039, 0x00599d16, 0x005acb94, 
  0x005bfbb8, 0x005d2d82, 0x005e60f5, 0x005f9613, 
  0x0060ccdf, 0x0062055b, 0x00633f89, 0x00647b6d, 
  0x0065b907, 0x0066f85b, 0x0068396a, 0x00697c38, 
  0x006ac0c7, 0x006c0719, 0x006d4f30, 0x006e9910, 
  0x006fe4ba, 0x00713231, 0x00728177, 0x0073d290, 
  0x0075257d, 0x00767a41, 0x0077d0df, 0x0079295a, 
  0x007a83b3, 0x007bdfed, 0x007d3e0c, 0x007e9e11
};
#else
#error invalid table size
#endif

inline unsigned int mask(int x)
{
	return (1U << x) - 1;
}

typedef union {
	float f;
	unsigned int i;
} fi;

float fmath_exp(float x) 
{
	const int s = FMATH_EXP_TABLE_SIZE;
	const int n = 1 << s; 
	const float a0 = n / logf(2.0);
	const float b0 = logf(2.0) / n;
	const float magic = (1 << 23) + (1 << 22); // to round

	float t = x * a0;
	t += magic;
	fi fi;
	fi.f = t;
	t = x - (t - magic) * b0;
	int u = ((fi.i + (127 << s)) >> s) << 23;
	unsigned int v = fi.i & mask(s);
	fi.i = u | kFmathExpTable[v];
	return (1.0f + t) * fi.f;
}

void fmath_exp4(float* RESTRICT y, const float* RESTRICT x) 
{
	const int s = FMATH_EXP_TABLE_SIZE;
	const int n = 1 << s; 
	const float a0 = n / logf(2.0);
	const float b0 = logf(2.0) / n;
	const float magic = (1 << 23) + (1 << 22); // to round

	float t0 = x[0] * a0;
	float t1 = x[1] * a0;
	float t2 = x[2] * a0;
	float t3 = x[3] * a0;
	t0 += magic;
	t1 += magic;
	t2 += magic;
	t3 += magic;
	fi fi0, fi1, fi2, fi3;
	fi0.f = t0;
	fi1.f = t1;
	fi2.f = t2;
	fi3.f = t3;
	t0 = x[0] - (t0 - magic) * b0;
	t1 = x[1] - (t1 - magic) * b0;
	t2 = x[2] - (t2 - magic) * b0;
	t3 = x[3] - (t3 - magic) * b0;
	int u0 = ((fi0.i + (127 << s)) >> s) << 23;
	int u1 = ((fi1.i + (127 << s)) >> s) << 23;
	int u2 = ((fi2.i + (127 << s)) >> s) << 23;
	int u3 = ((fi3.i + (127 << s)) >> s) << 23;
	unsigned int v0 = fi0.i & mask(s);
	unsigned int v1 = fi1.i & mask(s);
	unsigned int v2 = fi2.i & mask(s);
	unsigned int v3 = fi3.i & mask(s);
	fi0.i = u0 | kFmathExpTable[v0];
	fi1.i = u1 | kFmathExpTable[v1];
	fi2.i = u2 | kFmathExpTable[v2];
	fi3.i = u3 | kFmathExpTable[v3];
	y[0] = (1.0f + t0) * fi0.f;
	y[1] = (1.0f + t1) * fi1.f;
	y[2] = (1.0f + t2) * fi2.f;
	y[3] = (1.0f + t3) * fi3.f;
}

void fmath_exp8(float* y, const float* x) 
{
	const int s = FMATH_EXP_TABLE_SIZE;
	const int n = 1 << s; 
	const float a0 = n / logf(2.0);
	const float b0 = logf(2.0) / n;
	const float magic = (1 << 23) + (1 << 22); // to round

	float t0 = x[0] * a0;
	float t1 = x[1] * a0;
	float t2 = x[2] * a0;
	float t3 = x[3] * a0;
	float t4 = x[4] * a0;
	float t5 = x[5] * a0;
	float t6 = x[6] * a0;
	float t7 = x[7] * a0;
	t0 += magic;
	t1 += magic;
	t2 += magic;
	t3 += magic;
	t4 += magic;
	t5 += magic;
	t6 += magic;
	t7 += magic;
	fi fi0, fi1, fi2, fi3, fi4, fi5, fi6, fi7;
	fi0.f = t0;
	fi1.f = t1;
	fi2.f = t2;
	fi3.f = t3;
	fi4.f = t4;
	fi5.f = t5;
	fi6.f = t6;
	fi7.f = t7;
	t0 = x[0] - (t0 - magic) * b0;
	t1 = x[1] - (t1 - magic) * b0;
	t2 = x[2] - (t2 - magic) * b0;
	t3 = x[3] - (t3 - magic) * b0;
	t4 = x[3] - (t4 - magic) * b0;
	t5 = x[5] - (t5 - magic) * b0;
	t6 = x[6] - (t6 - magic) * b0;
	t7 = x[7] - (t7 - magic) * b0;
	int u0 = ((fi0.i + (127 << s)) >> s) << 23;
	int u1 = ((fi1.i + (127 << s)) >> s) << 23;
	int u2 = ((fi2.i + (127 << s)) >> s) << 23;
	int u3 = ((fi3.i + (127 << s)) >> s) << 23;
	int u4 = ((fi4.i + (127 << s)) >> s) << 23;
	int u5 = ((fi5.i + (127 << s)) >> s) << 23;
	int u6 = ((fi6.i + (127 << s)) >> s) << 23;
	int u7 = ((fi7.i + (127 << s)) >> s) << 23;
	unsigned int v0 = fi0.i & mask(s);
	unsigned int v1 = fi1.i & mask(s);
	unsigned int v2 = fi2.i & mask(s);
	unsigned int v3 = fi3.i & mask(s);
	unsigned int v4 = fi3.i & mask(s);
	unsigned int v5 = fi5.i & mask(s);
	unsigned int v6 = fi6.i & mask(s);
	unsigned int v7 = fi7.i & mask(s);
	fi0.i = u0 | kFmathExpTable[v0];
	fi1.i = u1 | kFmathExpTable[v1];
	fi2.i = u2 | kFmathExpTable[v2];
	fi3.i = u3 | kFmathExpTable[v3];
	fi4.i = u4 | kFmathExpTable[v4];
	fi5.i = u5 | kFmathExpTable[v5];
	fi6.i = u6 | kFmathExpTable[v6];
	fi7.i = u7 | kFmathExpTable[v7];
	y[0] = (1.0f + t0) * fi0.f;
	y[1] = (1.0f + t1) * fi1.f;
	y[2] = (1.0f + t2) * fi2.f;
	y[3] = (1.0f + t3) * fi3.f;
	y[4] = (1.0f + t4) * fi4.f;
	y[5] = (1.0f + t5) * fi5.f;
	y[6] = (1.0f + t6) * fi6.f;
	y[7] = (1.0f + t7) * fi7.f;
}

// Based on http://gallium.inria.fr/blog/fast-vectorizable-math-approx/

/* Relative error bounded by 1e-5 for normalized outputs
   Returns invalid outputs for nan inputs */
float expapprox(float val)
{
	/* Workaround a lack of optimization in gcc */
	const float exp_cst1 = 2139095040.f;
	const float exp_cst2 = 0.f;

	union {
		int i;
		float f;
	} xu, xu2;
	float val2, val3, val4, b;
	int val4i;
	val2 = 12102203.1615614f * val + 1065353216.f;
#if FMATH_EXP_DISABLE_RANGE_CHECK
	val4 = val2;
#else
	val3 = val2 < exp_cst1 ? val2 : exp_cst1;
	val4 = val3 > exp_cst2 ? val3 : exp_cst2;
#endif
	val4i = (int)val4;
	xu.i = val4i & 0x7F800000;
	xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
	b = xu2.f;
	return xu.f *
	       (0.509964287281036376953125f +
		b * (0.3120158612728118896484375f +
		     b * (0.1666135489940643310546875f +
			  b * (-2.12528370320796966552734375e-3f +
			       b * 1.3534179888665676116943359375e-2f))));
}

void expapprox4(float *RESTRICT dst, const float *RESTRICT src)
{
	// Manual code expansion of exparrpox() x 4.

	/* Workaround a lack of optimization in gcc */
	const float exp_cst1 = 2139095040.f;
	const float exp_cst2 = 0.f;

	const float kCoeff[5] = {
	    0.509964287281036376953125f, 0.3120158612728118896484375f,
	    0.1666135489940643310546875f, -2.12528370320796966552734375e-3f,
	    1.3534179888665676116943359375e-2f};

	union {
		int i;
		float f;
	} xu_0, xu_1, xu_2, xu_3, xu2_0, xu2_1, xu2_2, xu2_3;

	float val2_0, val2_1, val2_2, val2_3;
	float val3_0, val3_1, val3_2, val3_3;
	float val4_0, val4_1, val4_2, val4_3;
	float b0, b1, b2, b3;
	int val4i_0, val4i_1, val4i_2, val4i_3;

	val2_0 = 12102203.1615614f * src[0] + 1065353216.f;
	val2_1 = 12102203.1615614f * src[1] + 1065353216.f;
	val2_2 = 12102203.1615614f * src[2] + 1065353216.f;
	val2_3 = 12102203.1615614f * src[3] + 1065353216.f;

#if FMATH_EXP_DISABLE_RANGE_CHECK
	val4_0 = val2_0;
	val4_1 = val2_1;
	val4_2 = val2_2;
	val4_3 = val2_3;
#else
	val3_0 = val2_0 < exp_cst1 ? val2_0 : exp_cst1;
	val3_1 = val2_1 < exp_cst1 ? val2_1 : exp_cst1;
	val3_2 = val2_2 < exp_cst1 ? val2_2 : exp_cst1;
	val3_3 = val2_3 < exp_cst1 ? val2_3 : exp_cst1;

	val4_0 = val3_0 > exp_cst2 ? val3_0 : exp_cst2;
	val4_1 = val3_1 > exp_cst2 ? val3_1 : exp_cst2;
	val4_2 = val3_2 > exp_cst2 ? val3_2 : exp_cst2;
	val4_3 = val3_3 > exp_cst2 ? val3_3 : exp_cst2;
#endif

	val4i_0 = (int)val4_0;
	val4i_1 = (int)val4_1;
	val4i_2 = (int)val4_2;
	val4i_3 = (int)val4_3;

	xu_0.i = val4i_0 & 0x7F800000;
	xu_1.i = val4i_1 & 0x7F800000;
	xu_2.i = val4i_2 & 0x7F800000;
	xu_3.i = val4i_3 & 0x7F800000;

	xu2_0.i = (val4i_0 & 0x7FFFFF) | 0x3F800000;
	xu2_1.i = (val4i_1 & 0x7FFFFF) | 0x3F800000;
	xu2_2.i = (val4i_2 & 0x7FFFFF) | 0x3F800000;
	xu2_3.i = (val4i_3 & 0x7FFFFF) | 0x3F800000;

	b0 = xu2_0.f;
	b1 = xu2_1.f;
	b2 = xu2_2.f;
	b3 = xu2_3.f;

	const float c0 = kCoeff[0];
	const float c1 = kCoeff[1];
	const float c2 = kCoeff[2];
	const float c3 = kCoeff[3];
	const float c4 = kCoeff[4];

	dst[0] = xu_0.f * (c0 + b0 * (c1 + b0 * (c2 + b0 * (c3 + b0 * c4))));
	dst[1] = xu_1.f * (c0 + b1 * (c1 + b1 * (c2 + b1 * (c3 + b1 * c4))));
	dst[2] = xu_2.f * (c0 + b2 * (c1 + b2 * (c2 + b2 * (c3 + b2 * c4))));
	dst[3] = xu_3.f * (c0 + b3 * (c1 + b3 * (c2 + b3 * (c3 + b3 * c4))));
}

//
// -------------------------------------------------------------------------------------
//

#if FMATH_EXP_TEST

#ifndef WAIT_MICROSECONDS
#define WAIT_MICROSECONDS (10000)
#endif

#define TEST_NUM	  (4)

#include <stdio.h>

void validateExp(float retDiff[3], float beginValue, float endValue)
{
	int n = WAIT_MICROSECONDS / 250 / TEST_NUM; // 250 = emprically found value.
	union {
		int i;
		float f;
	} bv, ev, it;
	bv.f = beginValue;
	ev.f = endValue;
	float step = (endValue - beginValue) / n;

	it.i = bv.i;
	int count = 0;
	volatile float minDiff = 0.0f;
	volatile float maxDiff = 0.0f;
	volatile float aveDiff = 0.0f;
	// for (; it.i < ev.i; it.i++) { // <-- e-gcc can't compile this loop
	// ...
	float f = beginValue;
	for (f = beginValue; f < endValue; f += step) {

		float ref = expf(f);
		float ret = expapprox(f);
		float diff = fabsf(ref - ret) / ref;

		if (count == 0) {
			minDiff = diff;
			maxDiff = diff;
		} else {
			minDiff = (minDiff > diff) ? diff : minDiff;
			maxDiff = (maxDiff < diff) ? diff : maxDiff;
		}
		aveDiff += diff;
		count++;
	}

	aveDiff /= (float)count;

	retDiff[0] = aveDiff;
	retDiff[1] = minDiff;
	retDiff[2] = maxDiff;
}

void validateExp4(float retDiff[3], float beginValue, float endValue)
{
	int n = WAIT_MICROSECONDS / 250 / TEST_NUM; // 250 = emprically found value.
	union {
		int i;
		float f;
	} bv, ev, it;
	bv.f = beginValue;
	ev.f = endValue;
	float step = (endValue - beginValue) / n;

	it.i = bv.i;
	int count = 0;
	volatile float minDiff = 0.0f;
	volatile float maxDiff = 0.0f;
	volatile float aveDiff = 0.0f;
	// for (; it.i < ev.i; it.i++) { // <-- e-gcc can't compile this loop
	// ...
	float f = beginValue;
	for (f = beginValue; f < endValue; f += 4.0f * step) {

		float ref[4];
		ref[0] = expf(f);
		ref[1] = expf(f + step);
		ref[2] = expf(f + 2 * step);
		ref[3] = expf(f + 3 * step);
		float ret[4];
		float src[4] = {f, f + step, f + 2 * step, f + 3 * step};
		expapprox4(ret, src);

		for (int k = 0; k < 4; k++) {
			float diff = fabs(ref[k] - ret[k]) / (float)ref[k];

			if (count == 0) {
				minDiff = diff;
				maxDiff = diff;
			} else {
				minDiff = (minDiff > diff) ? diff : minDiff;
				maxDiff = (maxDiff < diff) ? diff : maxDiff;
			}
			aveDiff += diff;
			count++;
		}
	}

	aveDiff /= (float)count;

	retDiff[0] = aveDiff;
	retDiff[1] = minDiff;
	retDiff[2] = maxDiff;
}

void validateFmathExp(float retDiff[3], float beginValue, float endValue)
{
	int n = WAIT_MICROSECONDS / 250 / TEST_NUM; // 250 = emprically found value.
	union {
		int i;
		float f;
	} bv, ev, it;
	bv.f = beginValue;
	ev.f = endValue;
	float step = (endValue - beginValue) / n;

	it.i = bv.i;
	int count = 0;
	volatile float minDiff = 0.0f;
	volatile float maxDiff = 0.0f;
	volatile float aveDiff = 0.0f;
	// for (; it.i < ev.i; it.i++) { // <-- e-gcc can't compile this loop
	// ...
	float f = beginValue;
	for (f = beginValue; f < endValue; f += step) {

		float ref = expf(f);
		float ret = fmath_exp(f);
		float diff = fabsf(ref - ret) / ref;

		if (count == 0) {
			minDiff = diff;
			maxDiff = diff;
		} else {
			minDiff = (minDiff > diff) ? diff : minDiff;
			maxDiff = (maxDiff < diff) ? diff : maxDiff;
		}
		aveDiff += diff;
		count++;
	}

	aveDiff /= (float)count;

	retDiff[0] = aveDiff;
	retDiff[1] = minDiff;
	retDiff[2] = maxDiff;
}

void validateFmathExp4(float retDiff[3], float beginValue, float endValue)
{
	int n = WAIT_MICROSECONDS / 250 / TEST_NUM; // 250 = emprically found value.
	union {
		int i;
		float f;
	} bv, ev, it;
	bv.f = beginValue;
	ev.f = endValue;
	float step = (endValue - beginValue) / n;

	it.i = bv.i;
	int count = 0;
	volatile float minDiff = 0.0f;
	volatile float maxDiff = 0.0f;
	volatile float aveDiff = 0.0f;
	// for (; it.i < ev.i; it.i++) { // <-- e-gcc can't compile this loop
	// ...
	float f = beginValue;
	for (f = beginValue; f < endValue; f += 4.0f * step) {

		float ref[4];
		ref[0] = expf(f);
		ref[1] = expf(f + step);
		ref[2] = expf(f + 2 * step);
		ref[3] = expf(f + 3 * step);
		float ret[4];
		float src[4] = {f, f + step, f + 2 * step, f + 3 * step};
		fmath_exp4(ret, src);

		for (int k = 0; k < 4; k++) {
			float diff = fabs(ref[k] - ret[k]) / (float)ref[k];

			if (count == 0) {
				minDiff = diff;
				maxDiff = diff;
			} else {
				minDiff = (minDiff > diff) ? diff : minDiff;
				maxDiff = (maxDiff < diff) ? diff : maxDiff;
			}
			aveDiff += diff;
			count++;
		}
	}

	aveDiff /= (float)count;

	retDiff[0] = aveDiff;
	retDiff[1] = minDiff;
	retDiff[2] = maxDiff;
}

char outbuf[4096] SECTION("shared_dram");
int main(void)
{
	e_coreid_t coreid;
	unsigned int i;
	unsigned int num;
	unsigned int time_p;
	unsigned int time_c;
	unsigned int time_compare;
	unsigned *mailbox;
	unsigned int temp;

	float volatile in_sin;
	float volatile in_cos;
	float volatile in_sqt;
	float volatile in_ceil;
	float volatile in_log;
	float volatile in_exp;
	float volatile in_exp1;
	float volatile in_exp2;
	float volatile in_exp3;
	float volatile in_exp4;
	float volatile in_exp5;
	float volatile in_exp6;
	float volatile in_exp7;
	float re_f0, re_f1, re_f2, re_f3, re_f4, re_f5;

	in_exp = 1.88f;
	in_exp1 = 2.88f;
	in_exp2 = 3.88f;
	in_exp3 = 4.88f;
	in_exp4 = 4.88f;
	in_exp5 = 5.88f;
	in_exp6 = 6.88f;
	in_exp7 = 7.88f;
	mailbox = (unsigned *)0x6000;
	mailbox[0] = 0;
	mailbox[1] = 0xFFFFFFFF;
	mailbox[2] = 0xFFFFFFFF;
	mailbox[3] = 0xFFFFFFFF;

	// Who am I? Query the CoreID from hardware.
	coreid = e_get_coreid();
	sprintf(outbuf, "");

	const float in_exp_arr[8] = {in_exp,  in_exp1, in_exp2, in_exp3,
				     in_exp4, in_exp5, in_exp6, in_exp7};
	float out_exp_arr[8];

	// Get time waste on functions
	e_ctimer_set(E_CTIMER_0, E_CTIMER_MAX);
	time_p = e_ctimer_start(E_CTIMER_0, E_CTIMER_CLK);
	time_c = e_ctimer_get(E_CTIMER_0);
	e_ctimer_stop(E_CTIMER_0);
	time_compare = time_p - time_c;

	// expf() reference
	volatile unsigned int exp_ref_clocks = 0;
	{
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		volatile float ret = expf(in_exp);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		exp_ref_clocks = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf), "\nThe clock cycle count for "
						 "\"expf()\" (reference) is "
						 "%d.\n",
			exp_ref_clocks);
	}

	// fmath
	if (1) {
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		//volatile float ret = expapprox(in_exp);
		volatile float ret = fmath_exp(in_exp);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		temp = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf),
			"\nThe clock cycle count for \"fmath_exp()\" is %d.\n",
			temp);
	}

	if (1) { // fmath_exp4
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		fmath_exp4(out_exp_arr, in_exp_arr);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		// prevent compiler dead code optimization
		volatile float ret = out_exp_arr[0] + out_exp_arr[1] +
				     out_exp_arr[2] + out_exp_arr[3];

		temp = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf), "\nThe clock cycle count for "
						 "\"fmath_exp4()\" is %d (/4 = "
						 "%d).\n",
			temp, temp / 4);
	}

	// expapprox
	if (1) {
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		volatile float ret = expapprox(in_exp);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		temp = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf),
			"\nThe clock cycle count for \"expapprox()\" is %d.\n",
			temp);
	}

	if (1) { // expapprox4
		e_ctimer_set(E_CTIMER_1, E_CTIMER_MAX);
		e_ctimer_start(E_CTIMER_1, E_CTIMER_CLK);
		time_p = e_ctimer_get(E_CTIMER_1);

		expapprox4(out_exp_arr, in_exp_arr);

		time_c = e_ctimer_get(E_CTIMER_1);
		e_ctimer_stop(E_CTIMER_1);

		// prevent compiler dead code optimization
		volatile float ret = out_exp_arr[0] + out_exp_arr[1] +
				     out_exp_arr[2] + out_exp_arr[3];

		temp = time_p - time_c - time_compare;

		sprintf(outbuf + strlen(outbuf), "\nThe clock cycle count for "
						 "\"fmath_exp4()\" is %d (/4 = "
						 "%d).\n",
			temp, temp / 4);
	}
	// Validation
	{
		float diffs[3];
		validateExp(diffs, -30.0f, 30.0f);
		//validateExp4(diffs, -3.0f, 3.0f);
		//validateFmathExp(diffs, -30.0f, 30.0f);
		//validateFmathExp4(diffs, -30.0f, 30.0f);

		mailbox[0] = 1;
		mailbox[1] = *((unsigned int *)&diffs[0]); // ave
		mailbox[2] = *((unsigned int *)&diffs[1]); // min
		mailbox[3] = *((unsigned int *)&diffs[2]); // max
	}

	return EXIT_SUCCESS;
}
#endif
