****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
****************************************************************************************************
LATEST VERSION OF THE FILE:                                   materials.duke.edu/AFLOW/README_PROTO.TXT
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
PROTOTYPES LIB1

STR#           Symmetry  Name          TYPE      PROTO       SG       conc
-----------------------------------------------------------------------------------------------
FCC
A1,FCC,1,2        fcc/cF4   Fm-3m         A1        Cu          #225
A2,BCC,58,59      bcc/cI2   Im-3m         A2        W           #229
A3,HCP,115,117    hcp/hcp2  P6_3/mmc      A3        Mg          #194
A4,DIA,301,302    cF8       Fd-3m         A4        C (DIA)     #227  books
A5,betaSn,305,306 tI4       I41/amd       A5        betaSn      #141  Navy Database
A6,In,303,304     tI2       I4/mmm        A6        Indium      #139  Navy Database
A7,As,307,308     hR6       R-3m          A7        alphaAs     #166  Kolmogorov
A8
A9
A10,Hg,349,350    hR1       R-3m          A10       alphaHg     #166  ICSD

****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
PROTOTYPES BINARY

STR#           Pearson   Symmetry      Name     TYPE   PROTO      SG      notes
-----------------------------------------------------------------------------------------------
FCC
1,2            fcc/cF4   Fm-3m         A1        Cu              #225                              G(pic1N1)
1.cub,2.cub    fcc/cF4   Fm-3m         A1        Cu              #225                              G(pic1N1)
3              tP2       P4/mmm        L1_0      AuCu            #123     AB  along [001]          G(pic1N2,1)
4              hR2       R-3m          L1_1      CuPt            #166     AB  along [111]          G(pic1N2,2)
5,6            tI6       I4/mmm        beta1/2   AB2             #139     AB2 along [100]z         G(pic1N3,2,"C11_b")
7,8            oI6       Immm          -         MoPt2           #71      AB2 along [011]          G(pic1N3,1)
9,10           hP3       P-3m1         alpha1/2  AB2             #164     AB2 along [111]z         G(pic1N3,3,"C6")
11,12          oS8       Cmmm          "L1_3"(unofficial) AB3    #65      "(A.5B.5)B along [111]"  M(CdPt3,PdPt3,Pd3Pt) G(pic1N4,2)
13,15          tP4       P4/mmm        Z1/3      AB3             #123     AB3 along [001]z         M(PtTc3)  G(pic1N4,8)
14             tP4       P4/nmm        Z2        A2B2            #129     A2B2 along [001]z        G(pic1N4,3,"B11")
16,18          mS8       C2/m          ?         AB3             #12      AB3 along [113]          G(pic1N4,12)
17             mS8       C2/m          W2        A2B2            #12      A2B2 along [311]z        G(pic1N4,7)
19,21          oP4       Pmmm          Y1/3      AB3             #47      AB3  along [011]z        G(pic1N4,9)
20             oP4       Pmmm          Y2        A2B2            #47      A2B2 along [011]z        M(AgZr) G(pic1N4,4)
22,24          tI8       I4/mmm        D0_22     Al3Ti           #139     AB3  along [201]         G(pic1N4,11)
23             tI8       I41/amd       CH "40"   NbP             #141     A2B2 along [201]         G(pic1N4,6,NbP)
25,26          cP4       Pm3(bar)m     L1_2      AuCu3           #221     "(A.5B.5)B  along [001]" G(pic1N4,1)
27,29          hR4       R-3m          V1/3      AB3             #166     AB3  along [111]z        G(pic1N4,10)
28             hR4       R-3m          V2        A2B2            #166     A2B2 along [111]z        G(pic1N4,5)
BCC
58,59          bcc/cI2   Im-3m         A2        W               #229
58.cub,59.cub  bcc/cI2   Im-3m         A2        W               #229
60             oS4       Cmmm          D^19_2h   gamma-IrV       #65      AB along [101]
61             cP2       Pm-3m         B2        CsCl            #221     AB along [001]
62,63          hP3       P-3m1         ?         AB2             #164     AB2 along [211]
64,65          oF12      Fmmm          ?         AB2             #69      AB2 along [011] WARNI FIX NpPd
66,67          tI6       I4/mmm        C11_b     MoSi2           #139     AB2 along [001] WARNI FIX NpPd
68,69          hR4       R-3m          ?         AB3             #166     UNDEF
70,72          oS8       Cmmm          ?         AB3             #65      UNDEF
71             oS8       Cmma          ?         A2B2            #67      UNDEF
73,75          mP4       P2/m          ?         AB3             #10      UNDEF
74             mP4       P21/m         ?         A2B2            #11      UNDEF
76,78          tP4       P4/mmm        ?         AB3             #123     AB3 along [001]z
77             tP4       P4/nmm        B11       gamma-CuTi      #129     A2B2 along [001]z
79,81          oI8       Immm          ?         AB3             #71      UNDEF  M(MoTi3, Mo3Ti, Nb3Tc,MoTi3, Mo3Ti, Nb3Tc)
80             oI8       Imma          ?         A2B2            #74      UNDEF  M(MoTi)
82,83          tP4       P4/mmm        L6_0      CuTi3           #123
84,86          cF16      Fm-3m         D0_3      AlFe3/BiF3      #225
85             cF16      Fd-3m         B32       NaTl            #227
HCP
115,117        hcp/hcp2  P6_3/mmc      A3        Mg              #194
116            hP2       P_6m2         B_h       WC              #187
118,121        oP4       Pmm2          ?         AB3             #25      UNDEF
119            oP4       Pmmn          -         CuTe            #59
120            oP4       Pmma          B19       AuCd            #51
123            oI8       C2/m          ?         A2B2            #12      UNDEF
122,124        oI8       Imm2          ?         AB3             #44      UNDEF
125,127        hP4       P-6m2         ?         AB3             #187     UNDEF
126            hP4       P-3m1         ?         A2B2            #164     UNDEF M(RhRu)
128,132        mS12      Cm            ?         AB5             #8       UNDEF
129,134        mS12      C2/m          ?         A4B2            #12      BiHf2-134       #12  Hf-Paper
130,137        mS12      C2/m          ?         A4B2            #12
131,135        mS12      C2/m          ?         A4B2            #12
133,140        mS12      Cm            ?         A2B4            #8       UNDEF
136            mS12      Cm            ?         A3B3            #8       UNDEF
138,139        mS12      Cm            ?         A3B3            #8       UNDEF
141,145        oS12      Amm2          ?         AB5             #38
142,147        oS12      Cmcm          ?         A4B2            #63      M(RhRu2)
143,150        oS12      Cmcm          ?         A4B2            #63      mistaken for C49, beware
*144,148       oS12      Cmcm          ?         A4B2            #63      NOT RUN
*146,153       oS12      Amm2          ?         A2B4            #38      NOT RUN
*149,151       oS12      Amm2          ?         A3B3            #38      NOT RUN
*152,157       oS12      Amm2          ?         A3B3            #38      NOT RUN
*154,160       hP6       P-62m         ?         AB5             #189     NOT RUN
*155,158       hP6       P-62m         ?         A2B4            #189     NOT RUN
*156,159       hP6       Cmcm          ?         A2B4            #63      NOT RUN
*161,167       hR6       R32           ?         AB5             #155     NOT RUN
*162,165       hR6       R32           ?         A2B4            #155     NOT RUN
*163,166       hR6       C2/c          ?         A2B4            #15      NOT RUN
*164           hR6       C2            ?         A3B3            #5       NOT RUN
*168,169       hP6       P-6m2         ?         AB5             #187     NOT RUN
*170,175       hP6       P-6m2         ?         A2B4            #187     NOT RUN
*172           mS8       P-6m2         ?         A3B3            #187     NOT RUN
*171,176       hP6       P-3m1         ?         A2B4            #164     NOT RUN
*173,177       hP6       P63/mmc       ?         A2B4            #194     NOT RUN
*174           hP6       P3m1          ?         A3B3            #156     NOT RUN
LAVES AND OMEGAS
178,179        hP12      P6_3/mm       C14       MgZn2           #194    A4B8  Hex Laves Str
180,181        oP16      Pnma          D0_11     Fe3C            #62
182,183        cF24      Fd-3m         C15       Cu2Mg           #227    Cub Laves Phase
184,185        cP8       Pm3(bar)n     A15       Cr3Si           #223
186,187        hP8       P6_3/mmc      D0_19     Ni3Sn           #194
188,189        oS12      Cmcm          C49       ZrSi2           #63
190,191        hP3       P3m1      omega phase with Z=1/4 !      #164
192,193        oS8       Cmcm          B33B_f    CrB             #63
194            oS8       Cmcm                    ITl ITl         #63     Pauling
195            cP8       P2_13         B20       FeSi            #198    Pauling
196,197        tP4       P4/nmm O2               CdTi            #129    Pauling
198            oP24      Pbcm                    InTh            #57     Pauling
199            tI16      I4/mcm        B37       SeTl            #140    Pauling  

201            cF8       Fm-3m         B1        NaCl            #225
202,203        tI10      I4/mmm        D1_3      BaAl4           #139    NEW
204,205        hP6       P6/mmm        D2_d      CaCu5           #191
206,207        tI26      I4/mmm        D2_b      ThMn12          #139
208,209        hP9       P-62m         C22       Fe2P            #189
210,211        oP12      Pnma          C37       Co2Si           #62     AL1851.cif
212,213        hR19                              Th2Zn17         #166
214,215        cI28      I-43d         D7_3      Th3P4           #220
216,217        hP3       P6/mmm        C32       AlB2            #191
218            cF8       F-43m         B3        ZnS             #216
219            hP4       P63/mc        B4        ZnS             #186
220,221        hP4       P63/mmc       B8_1      NiAs            #194
222,223        hP16      P63/mcm       D8_8      Mn5Si3          #193    128601.cif
224,225        hP38      P63/mmc                 Th2Ni17         #194
226,227        cP8       Pm-3n         A15       Cr3Si           #223
228,229        tP6       P4/nmm        C38       Cu2Sb           #129
230,231        tI12                    C16       CuAl2           #140
232,233        tI6       I4/mmm        C11_b     MoSi2           #139
234,235        tl12      I4/mcm        C16       Al2Cu           #140    130617.cif
236,237        hP142     P63/mmc                 Cd58Y13         #194   ICSD_102087
****NO 236,237        cl184     Im-3          Cd6Y     (Y12Cd68)        #204    130045.cif
****NO 238,239        hp24      P-6m2                   Co3V            #187    130340.cif
238,239        hp24      P63/mmc       Pu3Al                     #194    Pauling FILE = Co3V
***240,241               P42/mnm       D8_b      CrFe            #136    130307.cif ???? HARD
242,243        oP24      PmnmO2                  NbPd3           #59     130379.cif
244,245                  P63/mmc       D0_24     Ni3Ti           #194    AL5957.cif
246,247        hP6       P63/mmc                 CaIn2           #194    AL5296.cif
248,249        tI12      I41/amd       C_c       ThSi2           #141    AL2272.cif
***250,251     tI6       I4/mmm        C11_a     CaC2            #139    AL1909.cif UNRUNNABLE
252,253        cF24      F-43m         C15_b     AuBe5           #216    AL3669.cif
254,255        tI32      I4/mcm                  W5Si3           #140    AL2650.cif
               hP5       P63/mmc       D5_2      La2O3           #194    130882.cif FRAC OCCUP
256            oP8       Pnma          B27       BFe             #62     AL1744.cif
257,258        oP6       Pnnm          C18       FeS2            #58     marcasite AL1657.cif
259,260        hR5       R-3m          C33       Bi2Te3          #166    120146.cif
261,262        cF96      Fd-3m                   Ti2Ni           #227
263,264        tI14      I4/mmm                  Ti3Cu4          #139    120654.cif
265,266        tP10                              Ti2Cu3          #129    IMPOSSIBLE
267,268        oP6       Pnnm          C18       FeS2            #52     marcasite AL1657.cif
269,270        hP3       P-3m1         C6        CdI2            #164    AL4899.cif
271,272        oS16      Cmcm                    Cd3Y            #63     18260.cif
273,274        hP6       P63/mmc       B8_2      Ni2In           #194    28249.cif
273p,274p      hP5       P-6m2      Fe1.5Ge1 partial of Ni2In    #187    Pauling S525640
275,276        hP6       P63/m                   Ni2Si           #176    AL5234.cif
277,278        oP8       PmmnO1        D0_a      betaCu3Ti       #59     AL1729.cif
279,280        tI16      I4/mmm        D0_23     Al3Zr           #139    Navy Database
281,282        cP12      Pa-3          C2        FeS2            #205    chris pyrite
283,284        oI12      Imma                    GdSi2(1.4)      #74     chris z=(0.0362,0.2138,0.6250)
285,286        tI10      I4/m          D1_a      MoNi4           #87     Navy Database
287,288        tI6       I4/mmm                  CuZr2           #139    Pauling file
289,290        cP4       Pm-3m         D0_9      alphaReO3       #221    Navy Database
291,292        tP4       P4/nmm        B10       PbO/PbS         #129    Navy Database
301,302        cF8       Fd-3m         A4        C (DIA)         #227    books
303,304        tI2       I4/mmm        A6        Indium          #139    Navy Database
305,306        tI4       I41/amd       A5        betaSn          #141    Navy Database
307,308        hR6       R-3m          A7        alphaAs         #166    Kolmo
309,310        cF32      Fm-3m                   CuPt7           #225    (Ca7Ge) Navy Database
311,312        tI18      I4/mmm                  NbNi8           #139    (Pt8Ti) Navy Database
313,314        tI18      I4/mmm                  V4Zn5           #139    Navy Database
315,316        hP24      P6_3/mmc      C36       MgNi2           #194    Navy Database
317,318        oS8       Cmca          A11       alphaGa         #64     Navy Database
319,320        cI58      I(-4)3m       A12       alphaMn         #217    Navy Database
321,322        cP20      P4_132        A13       betaMn          #213    Navy Database
323,324        hP3       P3_21-D_3^6   A8        gammaSe         #152    Navy Database
325,326        hP4       P63/mmc       A9        graphite        #194    Navy Database
327,328        tI2       I4/mmm        Aa        alphaPa         #139    Navy Database
329,330        tP30      P4_2/mnm      Ab        betaU           #136    Navy Database/Pauling
331,332        oP8       I4/mmm        Ac        alphaNp         #62     Navy Database/Pauling
333,334        tP4       P4/nmm        Ad        betaNp          #129    Navy Database/Pauling
335,336        hP1       P6/mmm        Af hex    HgSn_{6-10}     #191    Navy Database
337,338        todo Ag T-50 boron
339,340        todo Ah simple cubic
341,342        todo Ai betaPo
343,344        todo Ak alphaSe
345,346        todo Al betaSe
347,348        todo A3' alphaLa
349,350        hR1       R-3m          A10       alphaHg         #166    ICSD   A10.A/B
351,352        todo A14 In2
353,354        todo A16 alphaS
355,356        todo A17 black P
357,358        todo A20 alphaU
359,360        hP18     P63/mcm                  Ga4Ti5          #193    Pauling
361,362        tP20     P42/mnm                  Al2Zr3          #136    Pauling
363,364        hP7      P6/mmm                   Al3Zr4          #191    Pauling
365,366        oF40     F2dd                     Al3Zr2          #43     Pauling
367,368        cF112    Fm-3c                    NaZn13          #226    Pauling
369,370        cP7      Pm-3m                    CaB6            #221    Pauling
371,372        hP5      P-3m1          D5_19     Al3Ni2          #164    Navy Database
373,374        cI40     I-43m          Ga4Ni(Ga3.62Ni0.97)       #217    Navy Database (vacancies I23,     #197)
375,376        oS16     Cmmm                     Ga3Pt5          #65     Navy Database
377,378        mS44     C12/m1                   Ga9Ni13         #12     Pauling
379,380        cI112    Ia(-3)d                  Ga4Ni3          #230    Navy Database
381,382        tI24     I4_1/amd                 Ga2Hf           #141    Navy Database
383,384        cF116    Fm-m                     C6Mn23          #225    ICSD
385,386        oP44     Pnma                     Cu8Hf3          #62     ICSD
387,388        oS68     Cmca                     Ni10Zr7         #64     ICSD_240191
389,390        tP3      P4/mmm                   Hg2Pt           #123    ICSD_104339
391,392        oP20     Pnma                     B4Mg            #62     Pauling
393,394        oI64     Imma                     B7Mg            #74     Pauling
395,396        oP16     Pbam                     Ge3Rh5          #55     Pauling
397,398        tP14     P4/mbm                   In4Ti3          #127    Pauling
399,400        tP34     P4/nmm.O2      V8Sb9     V7.49Sb9        #129    Pauling
401,402        tI160    I-42d                    Li3B14          #122    Pauling/Kolmogorov
403,404        tP20     Pa/mbm                   LiB3            #127    Pauling/Kolmogorov
405,406        hR4      R-3m                     LiB-MS1         #166    Pauling  
407,408        hP8      P63/mmc                  LiB-MS2         #194    Pauling    

411,412        cI58     I-43m                    Re24Ti5         #217    Pauling
413,414        mS36     C12/m1                   Ni7Zr2          #12     Pauling
415,416        tI40     I4/m                     Pt11Zr9         #87     Pauling
417,418        oP20     Pnma                     Au4Zr           #62     Pauling
419,420        oI142    Immm                     Hf54Os17        #71     Pauling
421,422        hP9      P6_222         C40       CrSi2           #180    Pauling
423.424.t0     hP68     P6/m                  Ag51Gd14.o3.t0     #175    Pauling
423.424.t1     hP68     P6/m                  Ag51Gd14.o3.t1     #175    Pauling
423.424.t2     hP68     P6/m                  Ag51Gd14.o3.t2     #175    Pauling
423,424.t3     hP68     P6/m                  Ag51Gd14.o2.t3     #175    Pauling
423,424.t4     hP68     P6/m                  Ag51Gd14.o2.t4     #175    Pauling
423,424.t5     hP68     P6/m                  Ag51Gd14.o2.t5     #175    Pauling
423,424.t6     hP68     P6/m                  Ag51Gd14.o2.t6     #175    Pauling
423,424.t7     hP68     P6/m                  Ag51Gd14.o1.t7     #175    Pauling
423,424.t8     hP68     P6/m                  Ag51Gd14.o1.t8     #175    Pauling
423,424.t9     hP68     P6/m                  Ag51Gd14.o1.t9     #175    Pauling
425,426        cF184    Fd-3m                    Zn22Zr          #227    Pauling
427,428        oP24     Pnnm                     As2Ti           #58     Pauling
429,430        tI84     I4/mmm                   Ge10Ho11        #139    Pauling
431,432        hP48     P6_122                   Ir3Zr5          #178    Pauling
433,434        oP6      Pnnm                     CaCl2           #58     Pauling
435,436        cP5      P-43m                    CFe4            #215    Pauling
437,438        mS28     C12/c1                   C2Mn5           #15     Pauling
439,440        oP40     Pnma                     C3Mn7           #62     Pauling
441,442        hP20     P6_3mc                   Fe3Th7          #186    Pauling
443,444        hR24     R-3cr                    F3Fe            #167    Pauling
445,446        hR24     P321                     F3Fe            #150    Pauling  
447,448        mP4      P12_1/m1                 NiTi            #11     Pauling
449,450        cP6      Pn-3m          C3        Ag2O            #224    ICSD_281041 Cuprite
451,452        hP9      P-3            Bb        etaAg2Zn        #147    Navy Database
453,454        mS16     C2/m           D0_15     AlCl3           #12     ICSD_39566
455,456        hR8      R32            D0_14     AlF3            #155    ICSD_30274
457,458        hR10     R-3c           D5_1      alphaAl2O3      #167    ICSD_75479 Corundum
459,460        cI26     Im-3                     Al12W           #204    ICSD_58207
461,462        oP12     P212121                  Ag2Se           #19     ICSD_15213 Naumannite
463,464        hR26     R-3                      AlPd            #148    ICSD_58112
465,466        cP16     Pa-3           SC16      AlSb            #205    ICSD_41994
467,468        cI52     I-43m         gammaBrass Cu5Zn8          #217    Pauling
469,470        hR36     R-3m h                   BaPb3           #166    Pauling
471,472        hP8      Cmcm                     Hf3Sc*-h321     #63     Hf-Paper
473,474        hP6      P-62m                    Hf5Sc-h51       #189    Hf-Paper
475,476        mS12     C2/m                     BiHf2-134       #12     Hf-Paper
477,478        tP6      P4/mmm                   Hf5Pb-f63       #123    Hf-Paper
479,480        oS12     Cmmm                     HfPd5           #65     Hf-Paper
481,482        hR276    R-3ch                    Re25Zr21        #167.2  Pauling
483,484        tP6      P4/mmm        B11_3                      #123    Mike Mehl bcc_AAABBB^[001]
485,486        tP6      P4mm          B11_3'                     #99     Mike Mehl bcc_AABBAB^[001]
487,488        tP6      P4/mmm        Z3                         #123    Mike Mehl fcc_AAABBB^[001]
489,490        tP6      P4mm          Z3'                        #99     Mike Mehl fcc_AABBAB^[001]
491,492        oP102    Pmn2_1        Al13Co4                    #31     Pauling - Grin J_Al_Comp_206_243_1994
493,494        ms102    C1m1          Al12.13Co4   eta           #8      Pauling - Hurd Acta Cryst 15 441 1962
495,496        hP28     P6_3/mmc      Al5Co2     D8_11           #194    Pauling
497,498        mp22     P12_1/c1      Al9Co2                     #14     Pauling
499,500        tP10     P-42_1c                  Pd4Se           #114    Pauling
501,502        cF1192   Fd(-3)m                  NaCd2           #227    Samson, Nature 195, 259 (1962)
503,504        hP10     P6_3/m                   Er3Ru2          #176    Pauling
505,506        oP40     Pnma                     Pt3Sr7          #62     Pauling
507,508        cF120    Fm-3m                    Ir4Sc11         #225    Pauling
509,510        oP276    Pnna                     Ru25Y44         #52     Pauling
511,512        cP252    P4_132                   RuZn6           #213    Pauling
513,514        mP4      P21/m         B19'       NbRu-Beta''     #11     Eur.Phys.J.B 72 575(2009)
515,516        cF96     Fd-3m                    delta(CdNi)     #227.2  Pauling
517,518        hP3      P-3m1                    Cd2Ce           #164    Pauling
519,520        hP3      P6/mmm                   Hg2U            #191    Pauling
521,522        hP9      P-62m                    InMg2           #189    Pauling
523,524        oS16     Cmcm                     Cd3Er           #63     Pauling
525,526        hP8      P6_3/mmc                 CdMg3           #194    Pauling
527,528        hP24     P6_3/mmc                 CeNi3           #194    Pauling
529,530        oP32     Pnma                     CoSc3           #62     Pauling
531,532        tP4      P4/mmm                   Pb3Sr           #123    Pauling
533,534        hR36     R-3mh                    PuNi3           #166    Pauling
535,536        oP16     Pnma                     YZn3            #62     Pauling
537,538        hR4      R-3m                     Co2Y2*          #166    hcp-paper definition
539,540        mS12     C2/m                     Sc2Zr*          #12     hcp-paper definition
541,542        oI8      Immm                     Mo3Ti*-81       #71     hcp-paper/monster
543,544        oI8      Imma                     MoTi*-80        #74     hcp-paper/monster
545,546        hR3      P-3m1                    ReTi2*-81       #164    hcp-paper
547,548        tI6      I4/mmm                   Hf2Tl*-6        #139    hf-paper
549,550        oF12     Fmmm                     Be2Zn*-65       #69     hcp-paper
551,552        oI8      Imm2                     Re3Ru*-124      #44     hcp-paper
553,554        hP24     P63cm          D0_21     Cu2.82P         #185    Pauling
555,556        oP108    Pnma                     Mg2Au           #62     Pauling
557,558        cF16     Fm-3m                    BiF3            #225    Pearson's    crystal data
559,560        cI58     I-43m                    Al12Mg17        #217    Pauling
561,562        oS160    Cmcm                     MgAu3-x         #63     Pauling
563,564        oS64     Cmcm                     MgAu3+x         #63     Pauling
565,566        oF48     Fddd                     Mg2Cu           #70     Pauling
567,568        hP18     P-62c                    Mg2Ga           #190    Pauling
569,570        oP24     Pbam                     MgGa2           #55     Pauling
571,572        oI28     Ibam                     Mg5Ga2          #72     Pauling
573,574        tI28     I4/mmm                   Mg2Ga5          #139    Pauling
575,576        tI32     I41/a                    MgGa            #88     Pauling
577,578        hR72/24  R32                      Mg3Hg           #155    Pauling
579,580        hR48/16  R-3m                     Mg3In           #166    Pauling
581,582        cF408    F-43m                    Mg44Rh7         #216    Pauling
583,584        cF12     Fm-3m                    CaF2            #225    Pauling
585,586        hP8      P63/mmc        D0_18     Al3Ir           #194    Pauling
587,588        cF116    Fm-3m                    Mn23Th6         #225    Pauling
589,590        cP39     Pm-3                     Mg2Zn11         #200    Pauling
591,592        mS110    C2/m                     Mg4Zn7          #12     Pauling
593,594        cP20     P4_132                   Mg3Ru2          #213    Pauling
595,596        oS12     Cmcm                     Au2V            #63     Pauling
597,598        hP94     P63/mmc                  Sr9Mg38         #194    Pauling
600.ABCDE      tP30     P4_2/mnm          sigma  CrFe            #136    Pauling
600.XXXXX where X can be A,B (a total of 32 poscars)
611,612        cP140    Pm-3                     Rh13Sc57        #200    Pauling
613,614        hR13     Imm2                     Fe7W6           #166    Pauling
615,616        oP24     Pbca                     AuSn2           #61     Pauling
617,618        mS12     C2./m                    Bi2Pd           #12     Pauling
619,620        mS28     C2./m                    Bi2Pd5          #12     Pauling
621,622        oP16     Pnma                     Bi3Ni           #62     Pauling
623,624        oP32     Pnma                     Bi3Y5           #62     Pauling
625,626        cI120    Ia-3d                    Bi4Rh           #230    Pauling
627,628        hR48     R-3m                     BiMn3           #166    Pauling
629,630        oP16     Pmma                     BiPd3           #51     Pauling
631,632        tI12     I4/mmm                   BiTi2           #139    Pauling
633,634        mP12     P2_{1}./c                CoSb2           #14     Pauling
635,636        hR78     R3m                      Cu7Hg6          #160    Pearson (Gamma-Brass Cu5Zn8 Modified)
637,638        oI12     Imma                     Hg2K            #74     Pauling
639,640        cI52     I-43m                    Ir2Zn11         #217    Pauling files (Gamma-Brass Cu5Zn8 Modified Occupancy
641,642        tI32     I-4            D0_e      Ni3P            #82     Pauling
643.AB/BA      aP32     P-1                      KP15            #2      Pauling
644.AB/BA      tI116    I4_1/a                   Pu28Zr          #88     Pauling
645.AB/BA      cF52     Fm-3m                    B12U            #225    Pauling
646.AB/BA      oS68     Cmcm                     TiZn16          #63     Pauling
647.AB/BA      tI140    I4/mcm                   Rh2Y3           #140    Pauling
648.AB/BA      tP32     P4/ncc                   Rh3Pu5          #130    Pauling
649.AB/BA      tI32     I-42m                    alpha-SV3       #121    Pauling
650.AB         oP12     Pmma                     IrTa            #51     Pauling
651.AB/BA      oP36     Pnma                     Ge4Sm5          #62     Pauling
652.AB/BA      oP12     Pnma                     Cl2Pb           #62     Pauling
653.AB/BA      mP16     P12_1                    NbPt3           #11     Pauling
654.AB/BA      cI10     Im-3m                    Hg4Pt           #229    Pauling
655.AB/BA      hR14     R-3                      Pd4Pu3          #148    Pauling
656.AB/BA      hR15     R-3                      Er3Ni2          #148    Pauling
657.AB/BA      oI6      Immm                     Pd2Ti           #71     Pauling
658.AB/BA      tP8      P4/mmm                   Pd5Ti3          #123    Pauling
659.AB/BA      oS20     Cmcm                     Pd3Ti2          #63     Pauling
660.AB/BA      tP14     P4/mbm                   Hg5Mn2          #127    Pauling
661.AB/BA      oS32     Cmcm                     Pu3Pd5          #63     Pauling
662.AB/BA      oS12     Cmmm                     Ga2Zr           #65     G. Ghosh and M.Asta, Acta Mat. 53, 3225-3252 (2005).
663.AB/BA      mP30     P21/c                    FeP4_mP30       #14     Pauling
664.AB/BA      mS40     C2/C                     FeP4_mS40       #15     Pauling
665.AB/BA      tP16     P4/2mnm                  In3Ir           #136    Pauling
666.AB/BA      mS44     C2/m                     Fe6Ge5          #12     Pauling
667.AB/BA      mS42     C2/m                     Fe3Ga4          #12     Pauling
668.AB/BA      oS28     Cmcm                     Al6Mn           #63     Pauling
669.AB         oP8      Pnma                     AsFe            #62     Pauling
670.AB/BA      hP6      P6/mmm      B35          CoSn            #191    Pauling
671.AB/BA      cI32     Im-3        D0_2         As3Co           #204    Pauling
672.AB/BA      cI32     Im-3m                                    #229    Mike - 16 atom scell fcc, A15B1
673.AB/BA      cP16     Pm-3m                                    #221    Mike - 16 atom scell bcc, A15B1
674.A/B        hP3      P6/mmm                                   #191    Mike - omega phase (C32) in pure form
675.A/B        cP8      Pm-3n                                    #223    Mike - A15 phase in pure form
676.AB/BA      hP22     P6/mmc        Fe6.5Ge4-ideal(Fe14Ge8)    #194    Pauling
676p.AB/BA     oS42     Amm2          Fe6.5Ge4-part(Fe13Ge8)     #38     Pauling
677.AB/BA      hR26     R3m                      Al8Cr5          #160.2  ICSD_606753
678.AB/BA      tI24     I4mm                     Co5Ge7          #107    ICSD_43386
679.AB/BA      cP5      Pm-3m         E2_1                       #221    NAVY
680.AB/BA      tP12     P4/mmm                   beta-NW2        #123    NAVY (Hagg, sort of)
681.AB/BA      hR3      R3m                      alpha-MoS2      #160    NAVY
682.AB/BA      hR13     R-3m          D8_5r                      #166    NAVY
683.AB/BA      cI50     Im-3m                    Na13Cl12v7      #229    NAVY
684.A/B        cP8      P2_13                    alpha-N2        #198    NAVY
685.A/B        cP8      Pa-3                     alpha-N2p       #205    NAVY
686.A/B        hP4      P6_3/mmc                 beta-N2         #194    NAVY
687.A/B        hR16     R-3c                     epsilon-N2      #167    NAVY
688.A/B        tP4      P4_2/mnm                 gamma-N2        #136    NAVY
689.AB/BA      hP6      P6_3/mmc                 WN2_194af       #194    NAVY
690.AB/BA      mP12     P2_1/c                   Ag2Te           #14     NAVY
691.AB/BA      cF12     Fm-3m         C1                         #225    NAVY
692.AB/BA      hR3      R-3m          C19                        #166    NAVY
693.AB/BA      tP6      P4_2/mnm      C4                         #136    NAVY
694.AB/BA      hP6      P6_3/mmc      C7                         #194    NAVY
695.AB/BA      hP3      P-6m2                    WN2_187fg       #187    NAVY
696.AB/BA      hP6      P6_3/mmc                 WN2_194de       #194    NAVY
697.AB/BA      tP10     P4_2/mmc                 Hagg_N4W6       #131    NAVY
698.AB/BA      oP16     Pnma                     MoO3            #62     NAVY
699.AB/BA      hP12     P6/mmm                   O3W             #191    NAVY
700.AB/BA      oP16     Pnma                     P3Tc            #62     NAVY
701.AB/BA      hP5      P-3m1                    WN4_164ad       #164    NAVY
702.AB/BA      hP20     P6_3/mmc                 B4W             #194    NAVY
703.AB/BA      oP40     Pbca                     P4Re            #61     NAVY
704.AB/BA      hP14     P6_3/mmc     B8_1v2                      #194    NAVY
705.AB/BA      oS28     Cmca         B8_1v2p                     #64     NAVY
706.AB/BA      oS14     Cmmm                     Na4Cl3v1        #65     NAVY
707.AB/BA      cP7      P-43m                    C4W3            #215    NAVY
708.AB/BA      cP7      Pm-3m                    S3U4            #221    NAVY
709.AB/BA      mS14     C2/m                     O5Nb2           #12     NAVY
710.AB/BA      hP13     P-3m1        B8_1v3                      #164    NAVY
711.AB/BA      cF52     Fm-3m                    Na7Cl6_v3       #225    NAVY
712.AB/BA      cF60     Fm-3m                    Na7Cl8v1        #225    NAVY
713.AB/BA      hP15     P-3m1        B8_1v1                      #164    NAVY
714.AB/BA      hP15     P-6m2        B8_1v1p                     #187    NAVY
715.AB/BA      tP2      P4/mmm                   Na1Cl1v2        #123    NAVY
716.AB/BA      tP6      P4/mmm                   Na3Cl3v2        #123    NAVY
717.AB/BA      hP14     P3m1         B8_1v2pp                    #156    NAVY
718.AB/BA      hR6      R3m          B3v2                        #160    NAVY
719.AB/BA      cF56     F-43m                    Na7Cl7v2        #216    NAVY
720.AB/BA      cF8      Fm-3m        B1                          #225    NAVY
721.AB/BA      hP6      P-62m                    NTa             #189    NAVY
722.AB/BA      cP6      Pm-3m                    NbO             #221    NAVY
723.AB/BA      hP16     P6_3mc                   NbS             #186    NAVY
724.AB/BA      tP3      P4/mmm                   Hagg_NW2        #123    NAVY
725.AB/BA      tI12     I4_1/amd                 Mo2N            #141    NAVY
726.AB/BA      hP9      P-31m                    NV2             #162    NAVY
727.AB/BA      hP9      P-3m1                    NW2_Pearson     #164    NAVY
728.A/B        cP8      Pm-3n        beta-W                      #223    NAVY
729.AB/BA      oP76     Pnma                     Ge8Mn11         #62     Pauling
730.AB/BA      hP48     P-3c1                    NaTe3           #165    ICSD_61355
731.AB/BA      oP48     Pbcn                     NaTe            #60     ICSD_61382
732.AB/BA      mP40     P121/c1                  Sb3Sr2          #60     ICSD_32033
733.AB/BA      mS16     C12/m1                   SrAs3           #12     ICSD_41831
734.AB/BA      tI12     I4/mmm                   La2Sb           #139    ICSD_10442
735.AB/BA      oP32     Pnma                     Ga5Tm3          #62     ICSD_104009


**** BORIDES PROJECT
1000.AB/BA     cP7      Pm-3m                    B6Si1           #221    ICSD_20240    binary_boride_ICSD_frequency=21
1001.AB/BA     hP3      P6/mmm                   B2Zr1           #191    ICSD_30327    binary_boride ICSD frequency=17
1002.AB/BA     cF52     Fm-3m                    B12Y1           #225    ICSD_23860    binary_boride_ICSD_frequency=16
1003.AB/BA     hP3      P6/mmm                   B2Mg1           #191    ICSD_26675    binary_boride_ICSD_frequency=13
1004.AB/BA     tP20     P4/mbm                   B4La1           #127    ICSD_2360     binary_boride_ICSD_frequency=12
1005.AB/BA     oS8      Cmcm                     B1Ni1           #63     ICSD_26937    binary_boride_ICSD_frequency=8
1006.AB/BA     tI12     I4/mcm                   B1Mo2           #140    ICSD_24278    binary_boride_ICSD_frequency=8
1007.AB        cF8      F-43m                    As1B1           #216    ICSD_43871    binary_boride_ICSD_frequency=5
1008.AB/BA     oI14     Immm                     B4Cr3           #71     ICSD_33704    binary_boride_ICSD_frequency=5
1009.AB        cF8      Fm-3m                    B1Pu1           #225    ICSD_43661    binary_boride_ICSD_frequency=4
1010.AB/BA     tP10     P4/mbm                   B2V3            #127    ICSD_88317    binary_boride_ICSD_frequency=4
1011.AB/BA     hP3      P6/mmm                   B2Zr1           #191    ICSD_169458   binary_boride ICSD frequency=3
1012.AB/BA     oS20     Cmcm                     B3V2            #63     ICSD_79258    binary_boride_ICSD_frequency=3
1013.AB/BA     tP20     P4/mbm                   B4Ce1           #127    ICSD_612719   binary_boride ICSD frequency=3
1014.AB/BA     cF116    Fm-3m                    B6Fe23          #225    ICSD_54786    binary_boride_ICSD_frequency=3
1015.AB        hP2      P-6m2                    B1Tc1           #187    ICSD_168895   binary_boride_ICSD_frequency=2
1016.AB/BA     hP4      P6_{3}/mmc               B1Pt1           #194    ICSD_24363    binary_boride_ICSD_frequency=2
1017.AB/BA     hP6      P6_{3}/mmc               B2Re1           #194    ICSD_23871    binary_boride_ICSD_frequency=2
1018.AB/BA     oP6      Pmmn                     B2Ru1           #59     ICSD_31871    binary_boride_ICSD_frequency=2
1019.AB/BA     mS14     C2/m                     B5W2            #12     ICSD_20326    binary_boride_ICSD_frequency=2
1020.AB/BA     oP8      Pnma                     B1Fe1           #62     ICSD_391329   binary_boride ICSD frequency=2
1021.AB/BA     tI16     I4_{1}/amd               B1Mo1           #141    ICSD_24280    binary_boride_ICSD_frequency=2
1022.AB/BA     oS16     Cmcm                     B1Re3           #63     ICSD_43662    binary_boride_ICSD_frequency=2
1023.AB/BA     hP10     P6_{3}/mmc               B3Ru2           #194    ICSD_23715    binary_boride_ICSD_frequency=2
1024.AB/BA     oP16     Pnma                     B1Ni3           #62     ICSD_24306    binary_boride_ICSD_frequency=2
1025.AB/BA     hP20     P6_{3}mc                 B3Ru7           #186    ICSD_44343    binary_boride_ICSD_frequency=2
1026.AB        cP2      Pm-3m                    B1Tc1           #221    ICSD_168900   binary_boride_ICSD_frequency=1
1027.AB/BA     cF12     Fm-3m                    B1Be2           #225    ICSD_20384    binary_boride_ICSD_frequency=1
1028.AB/BA     hP3      P6/mmm                   B2Nb1           #191    ICSD_614887   binary_boride ICSD frequency=1
1029.AB/BA     oS8      Cmcm                     B1Rh1           #63     ICSD_150732   binary_boride_ICSD_frequency=1
1030.AB/BA     hP4      P6_{3}/mmc               B1Tc1           #194    ICSD_168897   binary_boride_ICSD_frequency=1
1031.AB/BA     hP4      P6_{3}/mmc               B1Tc1           #194    ICSD_168898   binary_boride_ICSD_frequency=1
1032.AB/BA     hR4      R-3m                     B2Li2           #166    Kolmogorov-Curtarolo LiB    binary_boride_ICSD_frequency=1
1033.AB/BA     oI10     Immm                     B4Cr1           #71     ICSD_24353    binary_boride ICSD frequency=1
1034.AB/BA     mS10     C2/m                     B4Mn1           #12     ICSD_15079    binary_boride_ICSD_frequency=1
1035.AB/BA     tI12     I4/mcm                   B1Fe2           #140    ICSD_160788   binary_boride ICSD frequency=1
1036.AB/BA     tI12     I4/mcm                   B1Fe2           #140    ICSD_160789   binary_boride_ICSD_frequency=1
1037.AB/BA     tI12     I-42m                    B1Fe2           #121    ICSD_16809    binary_boride_ICSD_frequency=1
1038.AB/BA     oP6      Pnnm                     B1Pd2           #58     ICSD_10487    binary_boride_ICSD_frequency=1
1039.AB/BA     hP6      P6_{3}/mmc               B1Pt2           #194    ICSD_615207   binary_boride_ICSD_frequency=1
1040.AB/BA     tI12     I4/mcm                   B1Ti2           #140    ICSD_189385   binary_boride_ICSD_frequency=1
1041.AB/BA     mS12     C2/m                     B2Mo1           #12     ICSD_418398   binary_boride ICSD frequency=1
1042.AB/BA     cF24     Fd-3m                    B2Sc1           #227    ICSD_106903   binary_boride_ICSD_frequency=1
1043.AB/BA     hR6      R-3m                     B4Mo2           #166    ICSD_39554    binary_boride_ICSD_frequency=1
1044.AB/BA     oI14     Immm                     B4Mn3           #71     ICSD_614732   binary_boride ICSD frequency=1
1045.AB/BA     oI14     Immm                     B4Nb3           #71     ICSD_76631    binary_boride ICSD frequency=1
1046.AB/BA     mS14     C2/m                     B5W2            #12     ICSD_20326    binary_boride_ICSD_frequency=1
1047.AB/BA     oP8      Pnma                     B1Co1           #62     ICSD_612863   binary_boride_ICSD_frequency=1
1048.AB/BA     hP8      P6_{3}22                 B1Fe3           #182    ICSD_184958   binary_boride_ICSD_frequency=1
1049.AB/BA     oP8      Pnma                     B1Li1           #62     ICSD_153291   binary_boride_ICSD_frequency=1
1050.AB/BA     oP8      Pnma                     B1Mn1           #62     ICSD_76630    binary_boride ICSD frequency=1
1051.AB/BA     hP8      P6_{3}/mmc               B2Li2           #194    Kolmogorov-Curtarolo LiB    binary_boride_ICSD_frequency=1
1052.AB/BA     hR8      R-3m                     B3Mo1           #166    ICSD_167734   binary_boride_ICSD_frequency=1
1053.AB/BA     hP8      P6_{3}/mmc               B3Re1           #194    ICSD_24361    binary_boride_ICSD_frequency=1
1054.AB/BA     hP10     P6_{3}/mmc               B3Ru2           #194    ICSD_108082   binary_boride_ICSD_frequency=1
1055.AB/BA     oP10     Pnnm                     B4Cr1           #58     ICSD_186851   binary_boride_ICSD_frequency=1
1056.AB/BA     oP10     Pnnm                     B4Fe1           #58     ICSD_425311   binary_boride ICSD frequency=1
1057.AB/BA     hP10     P6_{3}/mmc               B4Mo1           #194    ICSD_167735   binary_boride_ICSD_frequency=1
1058.AB/BA     oS22     Cmmm                     B6Ta5           #65     ICSD_68538    binary_boride_ICSD_frequency=1
1059.AB/BA     cF48     F-43m                    B11Li1          #216    ICSD_164841   binary_boride_ICSD_frequency=1
1060.AB/BA     hR12     R3m                      B11Li1          #160    ICSD_164842   binary_boride_ICSD_frequency=1
1061.AB/BA     cF48     F-43m                    B11Li1          #216    ICSD_164843   binary_boride_ICSD_frequency=1
1062.AB/BA     tI24     I-4m2                    B11Li1          #119    ICSD_164844   binary_boride_ICSD_frequency=1
1063.AB/BA     tI24     I-4m2                    B11Li1          #119    ICSD_164845   binary_boride_ICSD_frequency=1
1064.AB/BA     oF48     Fddd                     B1Cr2           #70     ICSD_188474   binary_boride_ICSD_frequency=1
1065.AB/BA     oP12     Pnma                     B1Rh2           #62     ICSD_24699    binary_boride_ICSD_frequency=1
1066.AB/BA     oP12     Pnma                     B2Fe1           #62     ICSD_425310   binary_boride_ICSD_frequency=1
1067.AB/BA     hP12     P6_{3}/mmc               B2W1            #194    ICSD_23716    binary_boride_ICSD_frequency=1
1068.AB/BA     tI26     I4/mmm                   B12Sc1          #139    ICSD_615424   binary_boride_ICSD_frequency=1
1069.AB/BA     mS28     C2/c                     B2Pd5           #15     ICSD_43513    binary_boride_ICSD_frequency=1
1070.AB/BA     mS28     C2/c                     B3Ni4           #15     ICSD_150561   binary_boride ICSD frequency=1
1071.AB/BA     mS28     C2/c                     B3Ni4           #15     ICSD_24308    binary_boride_ICSD_frequency=1
1072.AB/BA     hP14     P6_{3}/mmc               B5Mo2           #194    ICSD_167733   binary_boride_ICSD_frequency=1
1073.AB/BA     hP14     P6_{3}/mmc               B5W2            #194    ICSD_615696   binary_boride_ICSD_frequency=1
1074.AB/BA     hR15     R-3m                     B12Si3          #166    ICSD_615435   binary_boride ICSD frequency=1
1075.AB/BA     oP16     Pnma                     B1Co3           #62     ICSD_44339    binary_boride_ICSD_frequency=2
1076.AB/BA     oP16     Pnma                     B1Pd3           #62     ICSD_43514    binary_boride_ICSD_frequency=1
1077.AB/BA     tP16     P4_{2}/mmc               B4Mn4           #131    ICSD_74985    binary_boride_ICSD_frequency=1
1078.AB/BA     hP18     R-3m                     B2Mo1           #166    ICSD_418397   binary_boride ICSD frequency=1
1079.AB/BA     hP18     P6_{3}/mmc               B4Rh5           #194    ICSD_86395    binary_boride_ICSD_frequency=1
1080.AB/BA     hP20     P6_{3}/mmc               B4W1            #194    ICSD_43193    binary_boride_ICSD_frequency=1
1081.AB/BA     tP20     P4/mbm                   B4W1            #127    ICSD_615683   binary_boride_ICSD_frequency=1
1082.AB/BA     hP20     P6_{3}cm                 B9Li1           #185    ICSD_164625   binary_boride_ICSD_frequency=1
1083.AB/BA     oS46     Cmmm                     B20Na3          #65     ICSD_407159   binary_boride_ICSD_frequency=1
1084.AB/BA     oP28     Pnma                     B3Ni4           #62     ICSD_24307    binary_boride_ICSD_frequency=1
1085.AB/BA     oI64     Imma                     B3Si1           #74     ICSD_412621   binary_boride_ICSD_frequency=1
1086.AB/BA     oP38     Pbam                     B8Ru11          #55     ICSD_43663    binary_boride_ICSD_frequency=1


**** CARBIDE PROJECT
001.AA         fcc/cF4  Fm-3m         A1         Cu              #225    from gndstate
001.BB         fcc/cF4  Fm-3m         A1         Cu              #225    from gndstate
058.AA         bcc/cI2  Im-3m         A2         W               #229    from gndstate
058.BB         bcc/cI2  Im-3m         A2         W               #229    from gndstate
115.AA         hcp/hcp  P6_3/mmc      A3         Mg              #194    from gndstate
115.BB         hcp/hcp  P6_3/mmc      A3         Mg              #194    from gndstate
201.AB         cF8      Fm-3m         B1         NaCl            #225    from gndstate
061.AB         cP2      Pm-3m         B2         CsCl            #221    from gndstate
218.AB         cF8      F-43m         B3         ZnS             #216    from gndstate
219.AB         hP4      P63/mc        B4         ZnS             #186    from gndstate
116.AB         hP2      P_6m2         B_h        WC              #187    from gndstate
301.AA         cF8      Fd-3m         A4         C (DIA)         #227    BOOK
301.BB         cF8      Fd-3m         A4         C (DIA)         #227    BOOK
325.AA         hP4      P63/mmc       A9         C (GRA)         #194    BOOK
325.BB         hP4      P63/mmc       A9         C (GRA)         #194    BOOK
1200.AB        cF8      Fm-3m                    C1Ti1           #225    ICSD_1546     binary_carbide_ICSD_frequency=30
1201.AB        cF8      F-43m                    C1Si1           #216    ICSD_24217    binary_carbide_ICSD_frequency=20
1202.AB/BA     tI6      I4/mmm                   C2U1            #139    ICSD_14330    binary_carbide_ICSD_frequency=18
1203.AB/BA     cI40     I-43d                    C3Pr2           #220    ICSD_2444     binary_carbide_ICSD_frequency=17
1204.AB        hP2      P-6m2                    C1W1            #187    ICSD_15406    binary_carbide_ICSD_frequency=12
1205.AB        cP2      Pm-3m                    C1Sc1           #221    ICSD_181048   binary_carbide_ICSD_frequency=11
1206.AB/BA     cI132    Im-3                     C60K6           #204    ICSD_66879    binary_carbide_ICSD_frequency=5
1207.AB/BA     hR3      R-3m                     C1Y2            #166    ICSD_22283    binary_carbide_ICSD_frequency=4
1208.AB/BA     mS12     C2/c                     C2Ca1           #15     ICSD_54184    binary_carbide_ICSD_frequency=4
1209.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_16593    binary_carbide_ICSD_frequency=4
1210.AB/BA     tP70     P4/mnc                   C4Sc3           #128    ICSD_71145    binary_carbide_ICSD_frequency=4
1211.AB/BA     tI6      I4/mmm                   Ba1C2           #139    ICSD_56160    binary_carbide_ICSD_frequency=3
1212.AB/BA     hP3      P-3m1                    C1Nb2           #164    ICSD_43670    binary_carbide_ICSD_frequency=3
1213.AB/BA     oP12     Pbcn                     C1Mo2           #60     ICSD_65701    binary_carbide_ICSD_frequency=3
1214.AB/BA     cP12     Pa-3                     C2Th1           #205    ICSD_77470    binary_carbide_ICSD_frequency=3
1215.AB/BA     oP16     Pnma                     C2Rb2           #62     ICSD_51529    binary_carbide_ICSD_frequency=3 
1216.AB/BA     oP18     Pbam                     C5Y4            #55     ICSD_84869    binary_carbide_ICSD_frequency=3
1217.AB/BA     tP68     P-42_{1}c                C19Sc15         #114    ICSD_42631    binary_carbide_ICSD_frequency=1
1218.AB/BA     cF12     Fm-3m                    Be2C1           #225    ICSD_41567    binary_carbide_ICSD_frequency=2
1219.AB/BA     hP3      P-6m2                    C2Os1           #187    ICSD_168279   binary_carbide_ICSD_frequency=2
1220.AB/BA     oI8      Immm                     C2Li2           #71     ICSD_25705    binary_carbide_ICSD_frequency=2
1221.AB/BA     hP6      P6_{3}/mmc               C2Os1           #194    ICSD_168280   binary_carbide_ICSD_frequency=2
1222.AB/BA     oP12     Pbcn                     C1V2            #60     ICSD_9965     binary_carbide_ICSD_frequency=2
1223.AB/BA     mP22     P2_{1}/c                 C7Ho4           #14     ICSD_83381    binary_carbide_ICSD_frequency=2
1224.AB/BA     cF116    Fm-3m                    C6Cr23          #225    ICSD_69536    binary_carbide_ICSD_frequency=2
1225.AB/BA     hR3      R-3m                     Ba1C2           #166    ICSD_186576   binary_carbide_ICSD_frequency=1
1226.AB/BA     hP3      P-6m2                    C1Mo2           #187    ICSD_43669    binary_carbide_ICSD_frequency=1
1227.AB/BA     hP3      P-3m1                    C1Sc2           #164    ICSD_280743   binary_carbide_ICSD_frequency=1
1228.AB/BA     tI6      I4/mmm                   C2Re1           #139    ICSD_184661   binary_carbide_ICSD_frequency=1
1229.AB/BA     tP3      P4/mmm                   C2Re1           #123    ICSD_184664   binary_carbide_ICSD_frequency=1
1230.AB        hP4      P6_{3}mc                 C1Nb1           #186    ICSD_189092   binary_carbide_ICSD_frequency=1
1231.AB/BA     hP4      P6_{3}/mmc               C1Nb1           #194    ICSD_189093   binary_carbide_ICSD_frequency=1
1232.AB/BA     hP4      P6_{3}/mmc               C1Os1           #194    ICSD_168278   binary_carbide_ICSD_frequency=1
1233.AB/BA     hR4      R-3m                     C1Ru1           #166    ICSD_188285   binary_carbide_ICSD_frequency=1
1234.AB/BA     hP4      P6_{3}/mmc               C1Sc1           #194    ICSD_189088   binary_carbide_ICSD_frequency=1
1235.AB        hP4      P6_{3}mc                 C1Si1           #186    ICSD_24261    binary_carbide_ICSD_frequency=1
1236.AB/BA     tI10     I4/m                     C1Fe4           #87     ICSD_187144   binary_carbide_ICSD_frequency=1
1237.AB/BA     cP5      P-43m                    C1Fe4           #215    ICSD_44729    binary_carbide_ICSD_frequency=1
1238.AB/BA     oI12     Imm2                     B1C5            #44     ICSD_166554   binary_carbide_ICSD_frequency=1
1239.AB/BA     hP6      P3m1                     B1C5            #156    ICSD_166555   binary_carbide_ICSD_frequency=1
1240.AB/BA     oP6      Pnnm                     C1Co2           #58     ICSD_617391   binary_carbide_ICSD_frequency=1
1241.AB/BA     hP3      P6/mmm                   C1Fe2           #191    ICSD_162103   binary_carbide_ICSD_frequency=1
1242.AB/BA     oP6      Pnnm                     C1Fe2           #58     ICSD_76826    binary_carbide_ICSD_frequency=2
1243.AB/BA     hP6      P6_{3}/mmc               C1Re2           #194    ICSD_180376   binary_carbide_ICSD_frequency=1
1244.AB/BA     mS12     C2/m                     C2Ca1           #12     ICSD_54185    binary_carbide_ICSD_frequency=1
1245.AB/BA     mS12     C2/m                     C2Ca1           #12     ICSD_54188    binary_carbide_ICSD_frequency=1
1246.AB/BA     mS12     C2/m                     C2Ir1           #12     ICSD_181488   binary_carbide_ICSD_frequency=1
1247.AB/BA     tI12     I4/mcm                   C2Ir1           #140    ICSD_181489   binary_carbide_ICSD_frequency=1
1248.AB/BA     tP6      P4_{2}/mnm               C2Mg1           #136    ICSD_88057    binary_carbide_ICSD_frequency=1
1249.AB/BA     hP6      P6_{3}/mmc               C2Re1           #194    ICSD_184660   binary_carbide_ICSD_frequency=1
1250.AB/BA     tP6      P4_{2}/mmc               C2Si1           #131    ICSD_187721   binary_carbide_ICSD_frequency=1
1251.AB/BA     mS12     C2/c                     C2Th1           #15     ICSD_202306   binary_carbide_ICSD_frequency=1
1252.AB/BA     mS12     C2/c                     C2Th1           #15     ICSD_26568    binary_carbide_ICSD_frequency=1
1253.AB/BA     tP6      P4_{2}/mmc               C2Th1           #131    ICSD_77469    binary_carbide_ICSD_frequency=1
1254.AB/BA     hR7      R3m                      Al4C3           #160    ICSD_14397    binary_carbide_ICSD_frequency=1
1255.AB/BA     cP7      P-43m                    As4C3           #215    ICSD_163831   binary_carbide_ICSD_frequency=1
1256.AB/BA     cP7      Pm-3m                    C3Nb4           #221    ICSD_42758    binary_carbide_ICSD_frequency=1
1257.AB/BA     tP8      P-4m2                    B1C7            #115    ICSD_181952   binary_carbide_ICSD_frequency=1
1258.AB/BA     hP8      P3m1                     B1C7            #156    ICSD_181953   binary_carbide_ICSD_frequency=1
1259.AB/BA     cP8      P-43m                    B1C7            #215    ICSD_181954   binary_carbide_ICSD_frequency=1
1260.AB/BA     oP8      Pmm2                     B1C7            #25     ICSD_181955   binary_carbide_ICSD_frequency=1
1261.AB/BA     hR8      R3m                      B1C7            #160    ICSD_181956   binary_carbide_ICSD_frequency=1
1262.AB/BA     hP8      P6_{3}22                 C1Fe3           #182    ICSD_163151   binary_carbide_ICSD_frequency=1
1263.AB/BA     hP8      P6_{3}22                 C1Fe3           #182    ICSD_42542    binary_carbide_ICSD_frequency=1
1264.AB/BA     hP8      P6_{3}/mmc               C1Mo1           #194    ICSD_618298   binary_carbide_ICSD_frequency=1
1265.AB/BA     hR8      R-3c                     C1Ni3           #167    ICSD_17005    binary_carbide_ICSD_frequency=1
1266.AB/BA     cP8      P2_{1}3                  C1Os1           #198    ICSD_168277   binary_carbide_ICSD_frequency=1
1267.AB/BA     hP8      P6_{3}/mmc               C1Re1           #194    ICSD_618702   binary_carbide_ICSD_frequency=1
1268.AB        hP8      P6_{3}mc                 C1Si1           #186    ICSD_24170    binary_carbide_ICSD_frequency=1
1269.AB/BA     hP9      P-31m                    C1W2            #162    ICSD_162419   binary_carbide_ICSD_frequency=1
1270.AB/BA     hP9      P6/mmm                   C8Cs1           #191    ICSD_74641    binary_carbide_ICSD_frequency=1
1271.AB/BA     mP10     P2_{1}/m                 C1Fe4           #11     ICSD_187143   binary_carbide_ICSD_frequency=1
1272.AB        hP10	    P3m1                     C1Si1           #156    ICSD_107204   binary_carbide_ICSD_frequency=1
1273.AB        hR10    	R3m                      C1Si1           #160    ICSD_24168    binary_carbide_ICSD_frequency=1
1274.AB/BA     oS20     Cmcm                     C2Cr3           #63     ICSD_57009    binary_carbide_ICSD_frequency=1
1275.AB/BA     oP10     Pnnm                     C3Mg2           #58     ICSD_71941    binary_carbide_ICSD_frequency=1
1276.AB/BA     mP10     P2_{1}/c                 C4Ir1           #14     ICSD_181493   binary_carbide_ICSD_frequency=1
1277.AB/BA     oS20     Ccce                     C4Ir1           #68     ICSD_181496   binary_carbide_ICSD_frequency=1
1278.AB/BA     tP10     P4/nbm                   C4Ir1           #125    ICSD_181497   binary_carbide_ICSD_frequency=1
1279.AB/BA     mS22     C2/m                     C5Nb6           #12     ICSD_20695    binary_carbide_ICSD_frequency=1
1280.AB/BA     tI12     I-4m2                    B1C5            #119    ICSD_166553   binary_carbide_ICSD_frequency=1
1281.AB/BA     aP12     P-1                      B1C5            #2      ICSD_166556   binary_carbide_ICSD_frequency=1
1282.AB/BA     oP12     Pnma                     C1Fe2           #62     ICSD_187138   binary_carbide_ICSD_frequency=1
1283.AB/BA     hP12     P6_{3}/mmc               C1Mo1           #194    ICSD_44987    binary_carbide_ICSD_frequency=1
1284.AB/BA     oP12     Pbcn                     C1Mo2           #60     ICSD_43322    binary_carbide_ICSD_frequency=1
1285.AB/BA     oP12     Pnma                     C1Nb2           #62     ICSD_31973    binary_carbide_ICSD_frequency=1
1286.AB/BA     oP12     Pmma                     C1Nb2           #51     ICSD_653912   binary_carbide_ICSD_frequency=1
1287.AB        hP12 	  P6_{3}mc                 C1Si1           #186    ICSD_15325    binary_carbide_ICSD_frequency=1
1288.AB/BA     cF48     Fd-3m                    C1Ti2           #227    ICSD_77473    binary_carbide_ICSD_frequency=1
1289.AB/BA     mP12     P2_{1}/c                 C2Ca1           #14     ICSD_24074    binary_carbide_ICSD_frequency=1
1290.AB/BA     cP12     Pn-3m                    C2Ca1           #224    ICSD_31092    binary_carbide_ICSD_frequency=1
1291.AB/BA     aP12     P-1                      C2Ca1           #2      ICSD_66663    binary_carbide_ICSD_frequency=1
1292.AB/BA     hP12     P-62m                    C2Cs2           #189    ICSD_51533    binary_carbide_ICSD_frequency=1
1293.AB/BA     cP12     P4_{3}32                 C2Na1           #212    ICSD_186177   binary_carbide_ICSD_frequency=1
1294.AB/BA     hR13     R-3m                     C5Ti8           #166    ICSD_20822    binary_carbide_ICSD_frequency=1
1295.AB        hR14     R3m                      C1Si1           #160    ICSD_42513    binary_carbide_ICSD_frequency=1
1296.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_423885   binary_carbide_ICSD_frequency=1
1297.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_43194    binary_carbide_ICSD_frequency=1
1298.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_89328    binary_carbide_ICSD_frequency=1
1299.AB/BA     mS28     C2/c                     C2Mn5           #15     ICSD_618246   binary_carbide_ICSD_frequency=1
1300.AB/BA     mS28     C2/c                     C2Mn5           #15     ICSD_69535    binary_carbide_ICSD_frequency=1
1301.AB/BA     cI28     I-43d                    C3Sc4           #220    ICSD_42760    binary_carbide_ICSD_frequency=1
1302.AB/BA     hP14     P6_{3}/mmc               C6Eu1           #194    ICSD_169041   binary_carbide_ICSD_frequency=1
1303.AB/BA     hP14     P6_{3}/mmc               C6Yb1           #194    ICSD_601565   binary_carbide_ICSD_frequency=1
1304.AB/BA     hR15     R-3m                     B13C2           #166    ICSD_446      binary_carbide_ICSD_frequency=1
1305.AB/BA     hR15     R-3m                     B4C1            #166    ICSD_29093    binary_carbide_ICSD_frequency=1
1306.AB/BA     oP16     Pnma                     C1Cr3           #62     ICSD_617486   binary_carbide_ICSD_frequency=1
1307.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_167667   binary_carbide_ICSD_frequency=1
1308.AB/BA     tI32     I-4                      C1Fe3           #82     ICSD_187142   binary_carbide_ICSD_frequency=1
1309.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_99010    binary_carbide_ICSD_frequency=1
1310.AB        hR16    	R3m                      C1Si1           #160    ICSD_18136    binary_carbide_ICSD_frequency=1
1311.AB        hP16     P6_{3}mc                 C1Si1           #186    ICSD_24630    binary_carbide_ICSD_frequency=1
1312.AB/BA     cF64     Fd-3m                    C1Ti1           #227    ICSD_618950   binary_carbide_ICSD_frequency=1
1313.AB/BA     oP16     Pnma                     C2Cs2           #62     ICSD_51531    binary_carbide_ICSD_frequency=1
1314.AB/BA     tI32     I4_{1}/acd               C2K2            #142    ICSD_36142    binary_carbide_ICSD_frequency=1
1315.AB/BA     tI32     I4_{1}/acd               C2Na2           #142    ICSD_28066    binary_carbide_ICSD_frequency=1
1316.AB/BA     oP16     Pnma                     C3Ir1           #62     ICSD_181492   binary_carbide_ICSD_frequency=1
1317.AB/BA     oP16     Pbam                     C3Ir5           #55     ICSD_181485   binary_carbide_ICSD_frequency=1
1318.AB        hR18     R3m                      C1Si1           #160    ICSD_24631    binary_carbide_ICSD_frequency=1
1319.AB/BA     oF72     Fddd                     C8K1            #70     ICSD_70020    binary_carbide_ICSD_frequency=1
1320.AB/BA     oF72     Fddd                     C8Rb1           #70     ICSD_200563   binary_carbide_ICSD_frequency=1
1321.AB        hP20     P3m1                     C1Si1           #156    ICSD_38150    binary_carbide_ICSD_frequency=1
1322.AB/BA     oP20     Pnma                     C2Cr3           #62     ICSD_15086    binary_carbide_ICSD_frequency=1
1323.AB/BA     oP20     Pnma                     C2Cr3           #62     ICSD_653958   binary_carbide_ICSD_frequency=1
1324.AB/BA     hP20     P6_{3}mc                 C3Fe7           #186    ICSD_76830    binary_carbide_ICSD_frequency=1
1325.AB/BA     tP20     P-4c2                    C3Os2           #116    ICSD_168275   binary_carbide_ICSD_frequency=1
1326.AB/BA     oS40     P6_{3}/mmc               C4Ir1           #194    ICSD_181498   binary_carbide_ICSD_frequency=1
1327.AB/BA     cI40     Im-3m                    C7Ir3           #229    ICSD_181491   binary_carbide_ICSD_frequency=1
1328.AB/BA     mP22     P2_{1}/c                 C7Lu4           #14     ICSD_83382    binary_carbide_ICSD_frequency=1
1329.AB/BA     oP24     Pnma                     Ba1C2           #62     ICSD_168408   binary_carbide_ICSD_frequency=1
1330.AB/BA     cF112    Fd-3m                    C1Mn6           #227    ICSD_187037   binary_carbide_ICSD_frequency=1
1331.AB        hP28     P3m1                     C1Si1           #156    ICSD_28303    binary_carbide_ICSD_frequency=1
1332.AB/BA     cF116    Fm-3m                    C6Mn23          #225    ICSD_154720   binary_carbide_ICSD_frequency=1
1333.AB        hP32     P6_{3}mc                 C1Si1           #186    ICSD_42858    binary_carbide_ICSD_frequency=1
1334.AB/BA     tI64     I4_{1}/acd               C2K2            #142    ICSD_89528    binary_carbide_ICSD_frequency=1
1335.AB/BA     hP33     P3_{1}12                 C5V6            #151    ICSD_71098    binary_carbide_ICSD_frequency=1
1336.AB        hP36     P3m1                     C1Si1           #156    ICSD_42859    binary_carbide_ICSD_frequency=1
1337.AB/BA     tP36     P-4c2                    C5Ir4           #116    ICSD_181486   binary_carbide_ICSD_frequency=1
1338.AB        hP38     P3m1                     C1Si1           #156    ICSD_27635    binary_carbide_ICSD_frequency=1
1339.AB/BA     oP40     Pbcn                     C3Os2           #60     ICSD_168276   binary_carbide_ICSD_frequency=1
1340.AB        hP42     P3m1                     C1Si1           #156    ICSD_43827    binary_carbide_ICSD_frequency=1
1341.AB/BA     mP44     P2_{1}/c                 C7Ho4           #14     ICSD_74817    binary_carbide_ICSD_frequency=1
1342.AB/BA     cP60     P4_{3}32                 C0.88V1         #212    ICSD_60191    binary_carbide_ICSD_frequency=1
1343.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_22177    binary_carbide_ICSD_frequency=1
1344.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_43259    binary_carbide_ICSD_frequency=1
1345.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_85953    binary_carbide_ICSD_frequency=1
1346.AB/BA     tP68     P-42_{1}c                C19Er15         #114    ICSD_76817    binary_carbide_ICSD_frequency=3
1347.AB/BA     cP126    Pm-3n                    Ba3C60          #223    ICSD_70063    binary_carbide_ICSD_frequency=1
1348.AB/BA     cP126    Pm-3                     C60Sr3          #200    ICSD_75354    binary_carbide_ICSD_frequency=1


*** TERNARIES
T0001.A2BC/AB2C/ABC2              cF16  Fm(-3)m  AlCu2Mn-Heusler              #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C]
T0001[-tet[+,-]*%].A2BC/AB2C/ABC2 cF16  Fm(-3)m  AlCu2Mn-Heusler-tet*%        #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C] TET distorted along z +-*%
T0002.A2BC/AB2C/ABC2              cF16  F-43m    Cu1Li2Sn1-antiHeusler        #216    ICSD_15130 [Wyckoff A~C, 3 protos]
T0002[-tet[+,-]*%].A2BC/AB2C/ABC2 cF16  F-43m    Cu1Li2Sn1-antiHeusler-tet*%  #216    ICSD_15130 [Wyckoff A~C, 3 protos]
T0001.{ABC}                       cF16  Fm(-3)m  AlCu2Mn-Heusler              #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C ]
T0001[-tet[+,-]*%].{ABC}          cF16  Fm(-3)m  AlCu2Mn-Heusler-tet*%        #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C ] TET distorted along Z
T0002.{ABC}                       cF16  F-43m    Cu1Li2Sn1-antiHeusler        #216    ICSD_15130 [Wyckoff A~C, 3 protos]
T0002[-tet[+,-]*%].{ABC}          cF16  F-43m    Cu1Li2Sn1-antiHeusler-tet*%  #216    ICSD_15130 [Wyckoff A~C, 3 protos]  TET distorted along Z
T0003.{ABC}                       cF12  F-43m    AgAsMg-halfHeusler       C1b #216    ICSD_43819
T0003[-tet[+,-]*%].{ABC}          cF12  F-43m    AgAsMg-halfHeusler-tet*% C1b #216    ICSD_43819  TET distorted along Z
T0004.{ABC}    oP12     Pnma            Co1Ge1Mn1_ICSD_52968     #62     ICSD_52968
T0005.{ABC}    hP6      P63/mmc         Co1Ge1Mn1_ICSD_623495    #194    ICSD_623495
T0006.{ABC}    hP12     P63/mmc         Co3Ge1Mn2_ICSD_52972     #194    ICSD_52972
T0007.{ABC}    oP12     Pnma            Cl1La1Se1_ICSD_425686    #62     ICSD_425686
T0008.{ABC}    oP12     Pnma          Cd1Na1P1 Co2Si1_ICSD_44858 #62     ICSD_44858
T0009.{ABC}    cP5      Pm-3m                    perovskite      #221 
T0010.{ABC}    cF24     F-43m           Cu4MgSn       C15b       #216    ICSD_108141
T0011.{ABC}    tI10    I4/mmm           As2Ce1Pd2                #139    ICSD_604354

*** TERNARY SQS
TS001.ABC/BCA/CBA  oS64    Cmma                 SQS_L12          #67     SQS_L12

*** TERNARY FCC
TFCC001.ABC           hP3       P3m1         ternary_fcc_1        #156    3cell
TFCC002.ABC           oI6       Imm2         ternary_fcc_2        #44     3cell
TFCC003.ABC           tI6       I4mm         ternary_fcc_3        #107    3cell
TFCC004.ABC/BCA/CAB   hR4       R3m          ternary_fcc_4        #160    4cell
TFCC005.ABC/BCA/CAB   oP4       Pmm2         ternary_fcc_5        #25     4cell
TFCC006.ABC/BCA/CAB   mS8       Cm           ternary_fcc_6        #8      4cell
TFCC007.ABC/BCA/CAB   tP4       P4mm         ternary_fcc_7        #99     4cell
TFCC008.ABC/BCA/CAB   tI8       I-4m2        ternary_fcc_8        #119    4cell
TFCC009.ABC/BCA/CAB   hR4       R-3m         ternary_fcc_9        #166    4cell
TFCC010.ABC/BCA/CAB   oP4       Pmmm         ternary_fcc_10       #47     4cell
TFCC011.ABC/BCA/CAB   mS8       C2/m         ternary_fcc_11       #12     4cell
TFCC012.ABC/BCA/CAB   tP4       P4/mmm       ternary_fcc_12       #123    4cell
TFCC013.ABC/BCA/CAB   tI8       I4/mmm       ternary_fcc_13       #139    4cell
TFCC014.ABC/BCA/CAB   mP4       P2/m         ternary_fcc_14       #10     4cell
TFCC015.ABC/BCA/CAB   oS8       Cmmm         ternary_fcc_15       #65     4cell
TFCC016.ABC/BCA/CAB   tP4       P4/mmm       ternary_fcc_16       #123    4cell

*** TERNARY BCC
TBCC001.ABC           oF12      Fmm2         ternary_bcc_1        #42     3cell
TBCC002.ABC           hP3       P3m1         ternary_bcc_2        #156    3cell
TBCC003.ABC           tI6       I4mm         ternary_bcc_3        #107    3cell
TBCC004.ABC/BCA/CAB   oS8       Cmm2         ternary_bcc_4        #35     4cell
TBCC005.ABC/BCA/CAB   mP4       Pm           ternary_bcc_5        #6      4cell
TBCC006.ABC/BCA/CAB   tP4       P4mm         ternary_bcc_6        #99     4cell
TBCC007.ABC/BCA/CAB   oI8       Imm2         ternary_bcc_7        #44     4cell
TBCC008.ABC/BCA/CAB   cF16      F-43m        ternary_bcc_8        #216    4cell same as T0002 (anti-Heusler)
TBCC009.ABC/BCA/CAB   oS8       Cmmm         ternary_bcc_9        #65     4cell
TBCC010.ABC/BCA/CAB   mP4       P2/m         ternary_bcc_10       #10     4cell
TBCC011.ABC/BCA/CAB   tP4       P4/mmm       ternary_bcc_11       #123    4cell
TBCC012.ABC/BCA/CAB   oI8       Cmmm         ternary_bcc_12       #65     4cell
TBCC013.ABC/BCA/CAB   cF16      Fm-3m        ternary_bcc_13       #225    4cell same as T0001 (Heusler)
TBCC014.ABC/BCA/CAB   mS8       C2/m         ternary_bcc_14       #12     4cell
TBCC015.ABC/BCA/CAB   oP4       Pmmm         ternary_bcc_15       #47     4cell
TBCC016.ABC/BCA/CAB   tP4       P4/mmm       ternary_bcc_16       #123    4cell

*** QUATERNARY
Q0001.{ABCD}   cF40     Fm(-3)m         Al1F6K2Na1_ICSD_6027     #225    ICSD_6027 Elpasolite [Wyckoff all distinct, 24 protos]

*****************************************************************************************************************************************************
*****************************************************************************************************************************************************
ANRL PROTOTYPES
*****************************************************************************************************************************************************
ANRL - UNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
A_aP4_2_aci-001                                  aP4       P-1                 -         Cf                  #2        Cf (part 1)         
A_mP4_4_2a-001                                   mP4       P2_{1}              -         Te                  #4        High-Pressure Te (part 1)
A_mC12_5_3c-001                                  mC12      C2                  A19       Po                  #5        Po (part 1)         
A_mP8_10_2m2n-001                                mP8       P2/m                -         C                   #10       S-carbon (part 2)   
A_mP16_11_8e-001                                 mP16      P2_{1}/m            -         alpha-Pu            #11       alpha-Pu (part 1)    
A_mC16_12_4i-001                                 mC16      C2/m                -         C                   #12       M-carbon (part 2)   
A_mC34_12_ah3i2j-001                             mC34      C2/m                -         beta-Pu             #12       beta-Pu (part 1)     
A_mC4_12_i-001                                   mC4       C2/m                -         O2                  #12       alpha-O2 (part 1)    
A_mP84_13_21g-001                                mP84      P2/c                -         P                   #13       Monoclinic Phosphorus (part 1)
A_mP32_14_8e-001                                 mP32      P2_{1}/c            A_l       Se                  #14       betaSe (part 1)     
A_mP64_14_16e-001                                mP64      P2_{1}/c            A_k       Se                  #14       Se (part 1)         
A_mC24_15_2e2f-001                               mC24      C2/c                -         H                   #15       H-III (part 2)      
A_oP16_55_2g2h-001                               oP16      Pbam                -         C                   #55       R-carbon (part 2)   
A_oP8_62_2c-001                                  oP8       Pnma                A_c       alpha-Np            #62       alpha-Np (part 1)    
A_oC4_63_c-001                                   oC4       Cmcm                A20       alpha-U             #63       alpha-Uranium (part 1)
A_oC8_64_f-001                                   oC8       Cmce                A11       alpha-Ga            #64       alpha-Gallium (part 1)
A_oC8_64_f-002                                   oC8       Cmce                A17       P                   #64       Black Phosphorus (part 1)
A_oC8_64_f-003                                   oC8       Cmce                A14       I2                  #64       Molecular Iodine (part 1)
A_oF128_70_4h-001                                oF128     Fddd                A16       S                   #70       alpha-Sulfur (part 1)
A_oF8_70_a-001                                   oF8       Fddd                -         gamma-Pu            #70       gamma-Pu (part 1)   
A_tP12_96_ab-001                                 tP12      P4_{3}2_{1}2        -         Si                  #96       "ST12" of Si (part 1)
A_tP4_129_ac-001                                 tP4       P4/nmm              A_d       beta-Np             #129      beta-Np (part 1)     
A_tP50_134_b2m2n-001                             tP50      P4_{2}/nnm          A_g       B                   #134      T-50 Boron (part 1) 
A_tP30_136_bf2ij-001                             tP30      P4_{2}/mnm          A_b       beta-U              #136      beta-U (part 1)      
A_tP4_136_f-001                                  tP4       P4_{2}/mnm          -         N2                  #136      gamma-N2 (part 1)    
sigma_tP30_136_bf2ij-001                         tP30      P4_{2}/mnm          D8_b      CrFe                #136      sigma-CrFe (part 1)  
A_tP12_138_bi-001                                tP12      P4_{2}/ncm          -         C                   #138      C (part 2)          
A_tP16_138_j-001                                 tP16      P4_{2}/ncm          A18       Cl2                 #138      Cl2 (part 1)        
A_tI2_139_a-001                                  tI2       I4/mmm              A6        In                  #139      In (part 1)         
A_tI2_139_a-002                                  tI2       I4/mmm              A_a       Pa                  #139      alpha-Pa (part 1)    
A_tI4_139_e-001                                  tI4       I4/mmm              -         Si                  #139      Hypothetical BCT5 Si (part 1)
A_tI8_139_h-001                                  tI8       I4/mmm              -         C                   #139      Hypothetical Tetrahedrally Bonded Carbon with 4-Member Rings (part 1)
A_tI4_141_a-001                                  tI4       I4_{1}/amd          A5        Sn                  #141      betaSn (part 1)     
A_tI16_142_f-001                                 tI16      I4_{1}/acd          -         S                   #142      S-III (part 2)      
A_hP3_152_a-001                                  hP3       P3_{1}21            A8        Se                  #152      gammaSe (part 1)    
A_hP9_154_bc-001                                 hP9       P3_{2}21            -         S                   #154      S-II (part 2)       
A_hR105_166_bc9h4i-001                           hR105     R-3m                -         B                   #166      beta-Boron (part 1)  
A_hR12_166_2h-001                                hR12      R-3m                -         B                   #166      alpha-Boron (part 1) 
A_hR1_166_a-001                                  hR1       R-3m                A_i       Po                  #166      beta-Po (part 1)     
A_hR1_166_a-002                                  hR1       R-3m                A10       Hg                  #166      alpha-Hg (part 1)    
A_hR2_166_c-001                                  hR2       R-3m                A7        As                  #166      alpha-As (part 1)    
A_hR2_166_c-002                                  hR2       R-3m                -         C                   #166      Rhombohedral Graphite (part 1)
A_hR2_166_c-003                                  hR2       R-3m                -         O2                  #166      beta-O2 (part 1)     
A_hR3_166_ac-001                                 hR3       R-3m                C19       Sm                  #166      alpha-Sm (part 1)    
A_hP6_178_a-001                                  hP6       P6_{1}22            -         Sc                  #178      Sc-V (part 2)       
A_hP4_186_ab-001                                 hP4       P6_{3}mc            -         C                   #186      Buckled Graphite (part 1)
A_hP1_191_a-001                                  hP1       P6/mmm              A_f       gamma-HgSn6-10      #191      Simple Hexagonal Lattice (part 1)
A_hP2_194_c-001                                  hP2       P6_{3}/mmc          A3        Mg                  #194      Hexagonal Close Packed (part 1)
A_hP4_194_ac-001                                 hP4       P6_{3}/mmc          A3'       La                  #194      alpha-La (part 1)    
A_hP4_194_bc-001                                 hP4       P6_{3}/mmc          A9        C                   #194      Hexagonal Graphite structure (part 1)
A_hP4_194_f-001                                  hP4       P6_{3}/mmc          -         C                   #194      Lonsdaleite (part 1)
A_hP6_194_h-001                                  hP6       P6_{3}/mmc          -         C                   #194      Hypothetical Tetrahedrally Bonded Carbon with 3-Member Rings (part 1)
A_cP8_198_2a-001                                 cP8       P2_{1}3             -         N2                  #198      alpha-N2 (part 1)   
A_cF240_202_h2i-001                              cF240     Fm-3                -         C                   #202      FCC C60 Buckminsterfullerine (part 2)
A_cP240_205_10d-001                              cP240     Pa-3                -         C                   #205      Simple Cubic C60 Buckminsterfullerine (part 2)
A_cP8_205_c-001                                  cP8       Pa-3                -         N2                  #205      alpha-N2 (part 1)   
A_cI16_206_c-001                                 cI16      Ia-3                -         Si                  #206      BC8 (part 1)        
A_cP20_213_cd-001                                cP20      P4_{1}32            A13       Mn                  #213      beta-Mn (part 1)     
A_cI58_217_ac2g-001                              cI58      I-43m               A12       Mn                  #217      alpha-Mn (part 1)    
A_cI16_220_c-001                                 cI16      I-43d               -         Li                  #220      High Pressure cI16 Li (part 1)
A_cP1_221_a                                      cP1       Pm-3m               A_h       Po                  #221      alpha-Po (part 1)    
A_cP46_223_dik-001                               cP46      Pm-3n               -         Si                  #223      Si46 Clathrate (part 1)
A_cF4_225_a                                      cF4       Fm-3m               A1        Cu                  #225      Face-Centered Cubic (part 1), equivalent to 1,2,1.cub,2.cub [HTQC]
A_cF136_227_aeg-001                              cF136     Fd-3m               -         Si                  #227      Si34 Clathrate (part 1)
A_cF8_227_a                                      cF8       Fd-3m               A4        C                   #227      Diamond (part 1)    
A_cI2_229_a                                      cI2       Im-3m               A2        W                   #229      Body-Centered Cubic (part 1)

*****************************************************************************************************************************************************
ANRL - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
AB2_aP12_1_4a_8a-001.{AB}                        aP12      P1                  -         FeS2                #1        FeS2 (part 1)       
A2B_aP6_2_2i_i-001.{AB}                          aP6       P-1                 -         P2I4                #2        P2I4 (part 1)       
A2B_aP6_2_aei_i-001.{AB}                         aP6       P-1                 -         H2S                 #2        H2S (part 2)        
A2B_mP12_3_bc3e_2e-001.{AB}                      mP12      P2                  -         SiO2                #3        SiO2 (part 1)       
A8B5_mP13_6_a7b_3a2b-001.{AB}                    mP13      Pm                  -         Mo8P5               #6        Mo8P5 (part 2)      
AB_mP4_6_2b_2a-001.AB                            mP4       Pm                  -         FeNi                #6        FeNi (part 2)       
A2B_mP12_7_4a_2a-001.{AB}                        mP12      Pc                  -         H2S                 #7        H2S IV (part 2)     
A2B_mP18_7_6a_3a-001.{AB}                        mP18      Pc                  -         As2Ba               #7        As2Ba (part 2)      
A3B_mP16_7_6a_2a-001.{AB}                        mP16      Pc                  -         epsilon-WO3         #7        epsilon-WO3 (part 2)
A9B2_mP22_7_9a_2a-001.{AB}                       mP22      Pc                  -         Rh2Ga9              #7        Rh2Ga9 (part 2)     
A2B_mC144_9_24a_12a-001.{AB}                     mC144     Cc                  -         SiO2                #9        Monoclinic Low Tridymite (part 1)
A5B3_mC32_9_5a_3a-001.{AB}                       mC32      Cc                  -         alpha-P3N5          #9        alpha-P3N5 (part 2) 
AB3_mC16_9_a_3a-001.{AB}                         mC16      Cc                  -         H3Cl                #9        H3Cl (part 2)       
A2B_mP6_10_mn_bg-001.{AB}                        mP6       P2/m                -         delta-PdCl2         #10       delta-PdCl2 (part 2)
AB3_mP16_10_mn_3m3n-001.{AB}                     mP16      P2/m                -         H3Cl                #10       H3Cl (part 2)       
AB_mP6_10_en_am-001.AB                           mP6       P2/m                -         LiSn                #10       LiSn (part 2)       
AB_mP4_11_e_e-001.{AB}                           mP4       P2_{1}/m            -         NiTi                #11       NiTi2 (part 1)      
A5B2_mC14_12_a2i_i-001.{AB}                      mC14      C2/m                -         Au5Mn2              #12       Au5Mn2 (part 1)     
AB2_mC6_12_a_i-001.{AB}                          mC6       C2/m                C34       AuTe2               #12       Calaverite (part 1) 
AB3_mC16_12_g_ij-001.{AB}                        mC16      C2/m                D0_15     AlCl3               #12       AlCl3 (part 1)      
A2B_mP12_13_2g_ef-001.{AB}                       mP12      P2/c                -         H2S                 #13       H2S (part 2)        
A2B_mP12_14_2e_e-001.{AB}                        mP12      P2_{1}/c            C43       ZrO2                #14       Baddeleyite (part 1)
A2B_mP12_14_2e_e-002.{AB}                        mP12      P2_{1}/c            -         O2W                 #14       O2W binary oxide (R. Friedrich), ICSD #8217
A2B_mP12_14_2e_e-003.{AB}                        mP12      P2_{1}/c            -         O2V                 #14       O2V binary oxide (R. Friedrich), ICSD #647610
A2B_mP6_14_e_a-001.{AB}                          mP6       P2_{1}/c            -         gamma-PdCl2         #14       gamma-PdCl2 (part 2)
A7B8_mP120_14_14e_16e-001.{AB}                   mP120     P2_{1}/c            -         alpha-C7H8          #14       alpha-Toluene (part 2)
A2B5_mC28_15_f_e2f-001.{AB}                      mC28      C2/c                -         B2Pd5               #15       B2Pd5 (part 1)      
A2B5_mC28_15_f_e2f-002.{AB}                      mC28      C2/c                -         O5Sb2               #15       O5Sb2 binary oxide (R. Friedrich), ICSD #1422      
A2B_mC48_15_ae3f_2f-001.{AB}                     mC48      C2/c                -         SiO2                #15       Coesite (part 1)    
AB3_mC16_15_e_cf-001.{AB}                        mC16      C2/c                -         H3Cl                #15       H3Cl (part 2)       
AB_mC8_15_c_e-001.{AB}                           mC8       C2/c                B26       CuO                 #15       Tenorite (part 1)   
A2B_oP12_17_abe_e-001.{AB}                       oP12      P222_{1}            -         alpha-Ag2Se         #17       alpha-Naumannite (part 2)
AB3_oP16_18_ab_3c-001.{AB}                       oP16      P2_{1}2_{1}2        D0_{17}   BaS3                #18       BaS3 (part 1)
A2B_oP12_19_2a_a-001.{AB}                        oP12      P2_{1}2_{1}2_{1}    -         Ag2Se               #19       Naumannite (part 1) 
AB3_oP16_19_a_3a-001.{AB}                        oP16      P2_{1}2_{1}2_{1}    -         H3Cl                #19       H3Cl (part 2)       
A2B_oC24_20_abc_c-001.{AB}                       oC24      C222_{1}            -         SiO2                #20       Orthorhombic Tridymite (part 1)
AB2_oC6_21_a_k-001.{AB}                          oC6       C222                -         Ta2H                #21       Ta2H (part 2)       
AB_oF8_22_a_c-001.AB                             oF8       F222                -         FeS                 #22       FeS (part 2)        
A3B_oI32_23_ij2k_k-001.{AB}                      oI32      I222                -         H3S                 #23       H3S (part 2)        
AB_oP2_25_b_a-001.AB                             oP2       Pmm2                -         CdTe                #25       High-Pressure CdTe (part 1)
A2B_oP12_26_abc_ab-001.{AB}                      oP12      Pmc2_{1}            -         H2S                 #26       H2S (part 2)        
A2B_oP12_26_abc_ab-002.{AB}                      oP12      Pmc2_{1}            -         beta-SeO2           #26       beta-SeO2 (part 2)  
A5B_oP24_26_3a3b2c_ab-001.{AB}                   oP24      Pmc2_{1}            -         TlP5                #26       TlP5 (part 2)       
AB2_oP24_28_acd_2c3d-001.{AB}                    oP24      Pma2                C46       AuTe2               #28       Krennerite (part 1) 
A2B_oP12_29_2a_a-001.{AB}                        oP12      Pca2_{1}            -         ZrO2                #29       ZrO2 (part 2)       
AB2_oP12_29_a_2a-001.{AB}                        oP12      Pca2_{1}            -         FeS2                #29       Pyrite (part 2)     
A2B3_oP40_33_4a_6a-001.{AB}                      oP40      Pna2_{1}            -         kappa-Al2O3         #33       kappa-alumina (part 2)
AB_oP8_33_a_a-001.{AB}                           oP8       Pna2_{1}            -         CoAs                #33       Modderite (part 1)  
AB2_oP6_34_a_c-001.{AB}                          oP6       Pnn2                -         FeSb2               #34       FeSb2 (part 2)      
A2B_oC12_36_2a_a-001.{AB}                        oC12      Cmc2_{1}            C24       HgBr2               #36       HgBr2 (part 1)      
AB_oC8_36_a_a-001.{AB}                           oC8       Cmc2_{1}            -         HCl                 #36       HCl (part 2)        
A2B_oC12_38_de_ab-001.{AB}                       oC12      Amm2                -         Au2V                #38       Au2V (part 1)       
A2B3_oC40_39_2d_2c2d-001.{AB}                    oC40      Aem2                -         Ta3S2               #39       Ta3S2 (part 2)      
AB3_oC16_40_b_3b-001.{AB}                        oC16      Ama2                -         CeTe3               #40       CeTe3 (part 2)      
AB2_oC24_41_2a_2b-001.{AB}                       oC24      Aea2                C_e       PdSn2               #41       PdSn2 (part 1)      
AB4_oC20_41_a_2b-001.{AB}                        oC20      Aea2                D1_c      PtSn4               #41       PtSn4 (part 1)      
A10B3_oF52_42_2abce_ab-001.{AB}                  oF52      Fmm2                -         W3O10               #42       W3O10 (part 2)      
AB_oF8_42_a_a-001.AB                             oF8       Fmm2                -         BN                  #42       BN (part 2)         
AB2_oF72_43_ab_3b-001.{AB}                       oF72      Fdd2                C44       GeS2                #43       GeS2 (part 1)       
AB_oI4_44_a_b-001.AB                             oI4       Imm2                -         AsGa                #44       High-pressure GaAs (part 1)
A5B2_oP14_49_dehq_ab-001.{AB}                    oP14      Pccm                -         beta-Ta2O5          #49       beta-Ta2O5 (part 2) 
AB_oP4_51_e_f-001.AB                             oP4       Pmma                B19       AuCd                #51       beta'-AuCd (part 1) 
A2B_oP24_52_2e_cd-001.{AB}                       oP24      Pnna                -         GaCl2               #52       GaCl2 (part 2)      
A3B2_oP20_52_de_cd-001.{AB}                      oP20      Pnna                -         Sr2Bi3              #52       Sr2Bi3 (part 2)     
A2B_oP24_55_2g2h_gh-001.{AB}                     oP24      Pbam                -         GeAs2               #55       GeAs2 (part 2)      
A3B5_oP16_55_ch_agh-001.{AB}                     oP16      Pbam                -         Rh5Ge3              #55       Rh5Ge3 (part 2)     
A3B2_oP20_56_ce_e-001.{AB}                       oP20      Pccn                D5_11     Sb2O3               #56       Sb2O3 (part 1)      
AB_oP8_57_d_d-001.{AB}                           oP8       Pbcm                -         TlF                 #57       TlF-II, 57 (part 1) 
A2B_oP6_58_g_a-001.{AB}                          oP6       Pnnm                C50       alpha-PdCl2         #58       alpha-PdCl2 (part 2)
AB2_oP6_58_a_g-001.{AB}                          oP6       Pnnm                C35       CaCl2               #58       Hydrophilite (part 1)
AB2_oP6_58_a_g-002.{AB}                          oP6       Pnnm                -         etaFe2C             #58       eta-Fe2C (part 1)   
AB2_oP6_58_a_g-003.{AB}                          oP6       Pnnm                C18       FeS2                #58       Marcasite (part 1)  
A3B_oP8_59_bf_a-001.{AB}                         oP8       Pmmn                D0_a      TiCu3               #59       betaTiCu3 (part 1)  
AB_oP4_59_a_b-001.{AB}                           oP4       Pmmn                -         CuTe                #59       Vulcanite (part 1)  
A2B3_oP20_60_d_cd-001.{AB}                       oP20      Pbcn                -         Rh2S3               #60       Rh2S3 (part 2)      
A3B_oP32_60_3d_d-001.{AB}                        oP32      Pbcn                -         WO3                 #60       WO3 (part 2)        
A7B8_oP120_60_7d_8d-001.{AB}                     oP120     Pbcn                -         beta-C7H8           #60       beta-Toluene (part 2)
A2B_oP24_61_2c_c-001.{AB}                        oP24      Pbca                C21       TiO2                #61       Brookite (part 1)   
AB_oP16_61_c_c-001.{AB}                          oP16      Pbca                B_e       CdSb                #61       CdSb (part 1)       
AB_oP48_61_3c_3c-001.{AB}                        oP48      Pbca                -         Benzene             #61       Benzene (part 2)    
A2B3_oP20_62_2c_3c-001.{AB}                      oP20      Pnma                D5_{10}   Cr3C2               #62       Tongbaite (part 2)  
A2B_oP12_62_2c_c-001.{AB}                        oP12      Pnma                C37       Co2Si               #62       Co2Si (part 1)      
A2B_oP12_62_2c_c-002.{AB}                        oP12      Pnma                C25       HgCl2               #62       HgCl2 (part 1)      
A2B_oP12_62_2c_c-003.{AB}                        oP12      Pnma                C23       PbCl2               #62       Cotunnite (part 1)  
A2B_oP12_62_2c_c-004.{AB}                        oP12      Pnma                C29       SrH2                #62       SrH2 (part 2)      
A3B2_oP20_62_3c_2c-001.{AB}                      oP20      Pnma                D5_8      Sb2S3               #62       Stibnite (part 1)   
A3B7_oP40_62_cd_3c2d-001.{AB}                    oP40      Pnma                D10_1     C3Cr7               #62       C3Cr7 (part 1)      
A3B_oP16_62_cd_c-001.{AB}                        oP16      Pnma                D0_{20}   epsilon-NiAl3       #62       epsilon-NiAl3 (part 2)
A4B_oP20_62_2cd_c-001.{AB}                       oP20      Pnma                -         MgB4                #62       MgB4 (part 1)       
AB3_oP16_62_c_3c-001.{AB}                        oP16      Pnma                D0_{8}    MoO3                #62       Molybdite (part 2)  
AB3_oP16_62_c_cd-001.{AB}                        oP16      Pnma                D0_11     Fe3C                #62       Cementite (part 1)  
AB_oP8_62_c_c-001.{AB}                           oP8       Pnma                B16       GeS                 #62       GeS (part 1)        
AB_oP8_62_c_c-002.{AB}                           oP8       Pnma                B31       MnP                 #62       MnP (part 1)        
AB_oP8_62_c_c-003.{AB}                           oP8       Pnma                B27       FeB                 #62       FeB (part 1)        
AB_oP8_62_c_c-004.{AB}                           oP8       Pnma                B29       SnS                 #62       SnS (part 1)        
AB_oP8_62_c_c-005.{AB}                           oP8       Pnma                B14/B_{d} FeAs/NiSi           #62       Westerveldite/eta-NiSi (part 2/part 3); eta-NiSi is structurally similar to westerveldite
AB_oP8_62_c_c-006.{AB}                           oP8       Pnma                -         HgO                 #62       binary oxide (R. Friedrich), ICSD #14124 
A2B_oC12_63_2c_c-001.{AB}                        oC12      Cmcm                C49       ZrSi2               #63       ZrSi2 (part 1)      
A6B_oC28_63_efg_c-001.{AB}                       oC28      Cmcm                D2_{h}    MnAl6               #63       MnAl6 (part 2)      
AB_oC8_63_c_c-001.{AB}                           oC8       Cmcm                B33       CrB                 #63       CrB (part 1)        
AB_oC8_63_c_c-002.AB                             oC8       Cmcm                -         AB                  #63       sc derivative (Y. Lederer)
A2B_oC24_64_2f_f-001.{AB}                        oC24      Cmce                -         H2S                 #64       H2S (part 2)        
A3B5_oC16_65_ah_bej-001.{AB}                     oC16      Cmmm                -         Ga3Pt5              #65       Ga3Pt5 (part 1)     
AB3_oC8_65_a_bf-001.{AB}                         oC8       Cmmm                L1_3      CdPt3               #65       Predicted CdPt3 (part 1), fcc derivative (Y. Lederer), equivalent to 11,12 [HTQC]
AB_oC8_65_j_g-001.{AB}                           oC8       Cmmm                -         alpha-IrV           #65       alpha-IrV (part 1)   
A3B_oC64_66_gi2lm_2l-001.{AB}                    oC64      Cccm                -         H3S                 #66       H3S (part 2)        
A3B_oC64_66_kl2m_bdl-001.{AB}                    oC64      Cccm                -         beta-ThI3           #66       beta-ThI3 (part 2)  
AB_oC8_67_a_g-001.{AB}                           oC8       Cmme                -         alpha-FeSe          #67       alpha-FeSe (part 2) 
AB_oC8_67_a_g-002.{AB}                           oC8       Cmme                -         alpha-PbO           #67       alpha-PbO (part 2)  
AB4_oC20_68_a_i-001.{AB}                         oC20      Ccce                -         PdSn4               #68       PdSn4 (part 2)      
AB_oF8_69_a_b-001.AB                             oF8       Fmmm                B24       TlF                 #69       TlF (part 1)        
A2B_oF24_70_e_a-001.{AB}                         oF24      Fddd                C54       TiSi2               #70       TiSi2 (part 1)      
AB2_oF48_70_f_fg-001.{AB}                        oF48      Fddd                D1_{f}/C_{b} Mn2B/Mg2Cu       #70       Mn2B/Mg2Cu (part 2/part 3), two Strukturbericht prototypes are structurally similar       
A4B3_oI14_71_gh_cg-001.{AB}                      oI14      Immm                D7_{b}    Ta3B4               #71       Ta3B4 (part 2)      
AB2_oI6_71_a_g-001.{AB}                          oI6       Immm                -         MoPt2               #71       MoPt2 (part 1), fcc derivative (Y. Lederer), equivalent to 7,8 [HTQC]     
AB2_oI6_71_a_i-001.{AB}                          oI6       Immm                -         ReSi2               #71       ReSi2 (part 1)      
A2B_oI12_72_j_a-001.{AB}                         oI12      Ibam                C42       SiS2                #72       SiS2 (part 1)       
A2B_oI12_74_h_e-001.{AB}                         oI12      Imma                -         KHg2/CeCu2          #74       KHg2/CeCu2 (part 2/part 3) 
A4B_oI20_74_beh_e-001.{AB}                       oI20      Imma                D1_{b}    Al4U                #74       Al4U (part 2)       
A3B7_tP40_76_3a_7a-001.{AB}                      tP40      P4_{1}              -         Cs3P7               #76       Cs3P7 (part 2)      
A2B_tP48_77_8d_4d-001.{AB}                       tP48      P4_{2}              -         H2S                 #77       H2S III (part 2)    
AB2_tI48_80_2b_4b-001.{AB}                       tI48      I4_{1}              -         beta-NbO2           #80       beta-NbO2 (part 2)  
AB2_tP12_81_adg_2h-001.{AB}                      tP12      P-4                 -         GeSe2               #81       GeSe2 (part 2)      
A3B_tI32_82_3g_g-001.{AB}                        tI32      I-4                 D0_{e}    Ni3P                #82       Ni3P (part 2)       
A3B2_tP10_83_adk_j-001.{AB}                      tP10      P4/m                -         Ti2Ge3              #83       Ti2Ge3 (part 2)     
AB_tP16_84_cej_k-001.{AB}                        tP16      P4_{2}/m            B34       PdS                 #84       PdS (part 1)        
A2B_tP30_85_ab2g_cg-001.{AB}                     tP30      P4/n                -         SrBr2               #85       SrBr2 (part 2)      
AB3_tP32_86_g_3g-001.{AB}                        tP32      P4_{2}/n            -         Ti3P                #86       Ti3P (part 2)       
A4B5_tI18_87_h_ah-001.{AB}                       tI18      I4/m                -         Ti5Te4              #87       Ti5Te4 (part 1)     
AB4_tI10_87_a_h-001.{AB}                         tI10      I4/m                D1_a      Ni4Mo               #87       Ni4Mo (part 1)      
A4B_tI20_88_f_a-001.{AB}                         tI20      I4_{1}/a            -         ThCl4               #88       ThCl4 (part 2)      
AB2_tI96_88_2f_4f-001.{AB}                       tI96      I4_{1}/a            -         alpha-NbO2          #88       alpha-NbO2 (part 2) 
A2B_tP12_92_b_a-001.{AB}                         tP12      P4_{1}2_{1}2        -         SiO2                #92       alpha-Cristobalite (part 1)
A2B_tP12_92_b_a-002.{AB}                         tP12      P4_{1}2_{1}2        -         O2Ti                #92       O2Ti binary oxide (R. Friedrich), ICSD #161691, equivalent to alpha-O2Te (paratellurite, http://aflow.org/prototype-encyclopedia/A2B_tP12_92_b_a.TeO2.html, part 3)
A2B_tP36_96_3b_ab-001.{AB}                       tP36      P4_{3}2_{1}2        -         SiO2                #96       Keatite (part 1)    
A2B_tI12_98_f_a-001.{AB}                         tI12      I4_{1}22            -         CdAs2               #98       CdAs2 (part 2)      
A2B3_tP20_102_2c_b2c-001.{AB}                    tP20      P4_{2}nm            -         Gd3Al2              #102      Gd3Al2 (part 2)     
AB4_tP10_103_a_d-001.{AB}                        tP10      P4cc                -         NbTe4               #103      NbTe4 (part 2)      
A5B7_tI24_107_ac_abd-001.{AB}                    tI24      I4mm                -         Co5Ge7              #107      Co5Ge7 (part 2)     
AB_tI4_107_a_a-001.AB                            tI4       I4mm                -         GeP                 #107      GeP (part 2)        
A3B5_tI32_108_ac_a2c-001.{AB}                    tI32      I4cm                -         Sr5Si3              #108      Sr5Si3 (part 2)     
AB_tI8_109_a_a-001.AB                            tI8       I4_{1}md            -         NbAs                #109      NbAs (part 2)       
A2B_tP12_111_2n_adf-001.{AB}                     tP12      P-42m               -         MnF2                #111      MnF2 (part 2)       
AB_tP8_111_n_n-001.AB                            tP8       P-42m               -         VN                  #111      VN (part 2)
AB3_tP8_113_a_ce-001.{AB}                        tP8       P-42_{1}m           -         BaS3                #113      BaS3 (part 1)
A3B_tP32_114_3e_e-001.{AB}                       tP32      P-42_{1}c           -         SeO3                #114      SeO3 (part 2)       
A4B_tP10_114_e_a-001.{AB}                        tP10      P-42_{1}c           -         Pd4Se               #114      Pd4Se (part 2)      
A2B3_tP5_115_g_ag-001.{AB}                       tP5       P-4m2               -         Rh3P2               #115      Rh3P2 (part 2)      
AB2_tP12_115_j_egi-001.{AB}                      tP12      P-4m2               -         HgI2                #115      HgI2 (part 2)       
A2B3_tP20_116_bci_fj-001.{AB}                    tP20      P-4c2               -         Ru2Sn3              #116      Ru2Sn3 (part 2)     
A2B3_tP20_117_i_adgh-001.{AB}                    tP20      P-4b2               D5_{12}   beta-Bi2O3          #117      beta-Bi2O3 (part 2)
A3B_tP16_118_ei_f-001.{AB}                       tP16      P-4n2               -         RuIn3               #118      RuIn3 (part 2)      
A5B3_tP32_118_g2i_aceh-001.{AB}                  tP32      P-4n2               -         Ir3Ga5              #118      Ir3Ga5 (part 2)     
A3B_tI24_119_b2i_af-001.{AB}                     tI24      I-4m2               -         RbGa3               #119      RbGa3 (part 2)      
AB_tI4_119_c_a-001.AB                            tI4       I-4m2               -         GaSb                #119      GaSb (part 2)       
AB3_tP4_123_a_ce-001.{AB}                        tP4       P4/mmm              L6_0      CuTi3               #123      CuTi3 (part 1)      
AB_tP2_123_a_d-001.AB                            tP2       P4/mmm              L1_0      CuAu                #123      CuAu (part 1), fcc derivative (Y. Lederer), equivalent to 3 [HTQC] 
AB4_tP10_124_a_m-001.{AB}                        tP10      P4/mcc              -         NbTe4               #124      NbTe4 (part 2)      
A4B_tP10_125_m_a-001.{AB}                        tP10      P4/nbm              -         PtPb4               #125      PtPb4 (part 2)      
A2B3_tP10_127_g_ah-001.{AB}                      tP10      P4/mbm              D5_a      Si2U3               #127      Si2U3 (part 1)      
A4B_tP20_127_ehj_g-001.{AB}                      tP20      P4/mbm              D1_{e}    ThB4                #127      ThB4 (part 2)       
A2B_tP6_129_ac_c-001.{AB}                        tP6       P4/nmm              C38       Cu2Sb               #129      Cu2Sb (part 1)      
AB_tP4_129_a_c-001.{AB}                          tP4       P4/nmm              B10       PbO                 #129      PbO (part 1)        
AB_tP4_129_c_c-001.{AB}                          tP4       P4/nmm              B11       gamma-CuTi          #129      gamma-CuTi (part 1)  
AB_tP4_129_c_c-002.AB                            tP4       P4/nmm              -         AB                  #129      fcc derivative (Y. Lederer), equivalent to 14 [HTQC]
A5B3_tP32_130_cg_cf-001.{AB}                     tP32      P4/ncc              -         Ba5Si3              #130      Ba5Si3 (part 2)     
AB_tP4_131_c_e-001.{AB}                          tP4       P4_{2}/mmc          B17       PtS                 #131      PtS (part 1)        
AB_tP4_131_c_e-002.{AB}                          tP4       P4_{2}/mmc          -         PdO                 #131      PdO binary oxide (R. Friedrich), ICSD #26598        
AB3_tP32_133_h_i2j-001.{AB}                      tP32      P4_{2}/nbc          -         beta-V3S            #133      beta-V3S (part 2)   
A2B_tP24_135_gh_h-001.{AB}                       tP24      P4_{2}/mbc          C47       SeO2                #135      Downeyite (part 2)  
A2B_tP6_136_f_a-001.{AB}                         tP6       P4_{2}/mnm          C4        TiO2                #136      CrSi2 (part 1)      
AB_tP8_136_g_f-001.{AB}                          tP8       P4_{2}/mnm          -         beta-BeO            #136      betaBeO (part 1)    
A2B3_tP40_137_cdf_3g-001.{AB}                    tP40      P4_{2}/nmc          D5_{9}    Zn3P2               #137      Zn3P2 (part 2)      
A2B_tP6_137_d_a-001.{AB}                         tP6       P4_{2}/nmc          -         ZrO2                #137      ZrO2 (part 2)       
AB2_tP6_137_a_d-001.{AB}                         tP6       P4_{2}/nmc          C13       HgI2                #137      HgI2 (part 2)       
A12B_tI26_139_fij_a-001.{AB}                     tI26      I4/mmm              D2_b      Mn12Th              #139      Mn12Th (part 1)     
A2B_tI6_139_d_a-001.{AB}                         tI6       I4/mmm              L\'2      ThH2                #139      ThH2 (part 1)       
A3B_tI16_139_cde_e-001.{AB}                      tI16      I4/mmm              D0_23     Al3Zr               #139      Al3Zr (part 1)      
A3B_tI8_139_bd_a-001.{AB}                        tI8       I4/mmm              D0_22     Al3Ti               #139      Al3Ti (part 1)      
A3B_tI8_139_bd_a-002.{AB}                        tI8       I4/mmm              -         AB3                 #139      fcc derivative (Y. Lederer) equivalent to 22,24 [HTQC]
A4B5_tI18_139_i_ah-001.{AB}                      tI18      I4/mmm              -         V4Zn5               #139      V4Zn5 (part 1)      
A4B_tI10_139_de_a-001.{AB}                       tI10      I4/mmm              D1_3      Al4Ba               #139      Al4Ba (part 1)      
A8B_tI18_139_hi_a-001.{AB}                       tI18      I4/mmm              -         Pt8Ti               #139      Pt8Ti (part 1)      
AB2_tI6_139_a_e-001.{AB}                         tI6       I4/mmm              C11_b     MoSi2               #139      MoSi2 (part 1)      
AB2_tI6_139_a_e-002.{AB}                         tI6       I4/mmm              -         AB2                 #139      fcc derivative (Y. Lederer), equivalent to 5,6,6b [HTQC]
AB2_tI6_139_a_e-003.{AB}                         tI6       I4/mmm              -         BaO2                #139      BaO2 binary oxide (R. Friedrich), ICSD #24248 
AB2_tI6_139_a_e-004.{AB}                         tI6       I4/mmm              C11_{a}   KO2 (C2Ca-I)        #139      KO2 binary oxide (R. Friedrich), ICSD #38245, equivalent to 250,251,1202.AB/BA,1211.AB/BA [HTQC], equivalent to C2Ca-I (C11_{a}, http://aflow.org/prototype-encyclopedia/A2B_tI6_139_e_a.html, part 3)
AB_tI8_139_e_e-001.{AB}                          tI8       I4/mmm              D3_{1}    Hg2Cl2              #139      Calomel (part 2)    
AB_tI8_139_e_e-002.AB                            tI8       I4/mmm              -         AB                  #139      sc derivative (Y. Lederer)
A2B_tI12_140_h_a-001.{AB}                        tI12      I4/mcm              C16       Al2Cu               #140      Khatyrkite (part 1) 
A3B5_tI32_140_ah_bk-001.{AB}                     tI32      I4/mcm              D8_{m}    W5Si3               #140      W5Si3 (part 2)      
A3B5_tI32_140_ah_cl-001.{AB}                     tI32      I4/mcm              D8_{l}    Cr5B3               #140      Cr5B3 (part 2)      
AB3_tI16_140_b_ah-001.{AB}                       tI16      I4/mcm              D0_c      SiU3                #140      SiU3 (part 1)       
AB_tI16_140_ab_h-001.{AB}                        tI16      I4/mcm              B37       SeTl                #140      SeTl (part 1)       
A2B3_tI80_141_ceh_3h-001.{AB}                    tI80      I4_{1}/amd          -         In2S3               #141      betaIn2S3 (part 1)  
A2B_tI12_141_e_a-001.{AB}                        tI12      I4_{1}/amd          C5        TiO2                #141      Anatase (part 1)    
A2B_tI12_141_e_a-002.{AB}                        tI12      I4_{1}/amd          C_{c}     alpha-ThSi2         #141      alpha-ThSi2 (part 2)
A2B_tI24_141_2e_e-001.{AB}                       tI28      I4_{1}/amd          -         Ga2Hf               #141      Ga2Hf (part 1)      
A3B4_tI28_141_ad_h-001.{AB}                      tI28      I4_{1}/amd          -         Mn3O4               #141      Hausmannite (part 1)
AB_tI16_141_e_e-001.{AB}                         tI16      I4_{1}/amd          B_g       MoB                 #141      MoB (part 1)        
AB_tI8_141_a_b-001.AB                            tI8       I4_{1}/amd          "40"      NbP                 #141      NbP (part 1)        
AB_tI8_141_a_b-002.AB                            tI8       I4_{1}/amd          -         AB                  #141      fcc derivative (Y. Lederer), equivalent to 23 [HTQC]
AB2_hP12_143_cd_ab2d-001.{AB}                    hP12      P3                  -         MoS2                #143      MoS2 (part 2)       
A4B_hP15_144_4a_a-001.{AB}                       hP15      P3_{1}              -         IrGe4               #144      IrGe4 (part 2)      
AB_hP6_144_a_a-001.AB                            hP6       P3_{1}              -         ZnTe                #144      ZnTe (part 2)
A2B_hP9_147_g_ad-001.{AB}                        hP9       P-3                 B_b       AgZn                #147      zetaAgZn (part 1)   
A2B_hR18_148_2f_f-001.{AB}                       hR18      R-3                 -         beta-PdCl2          #148      beta-PdCl2 (part 2) 
AB3_hR8_148_c_f-001.{AB}                         hR8       R-3                 D0_5      BiI3                #148      BiI3 (part 1)       
AB_hR16_148_cf_cf-001.{AB}                       hR16      R-3                 -         C8H8                #148      Solid Cubane (part 1)
AB_hR26_148_b2f_a2f-001.AB                       hR26      R-3                 -         PdAl                #148      PdAl (part 1)       
AB3_hP24_149_acgi_3l-001.{AB}                    hP24      P312                -         Ti3O                #149      Ti3O (part 2)       
A2B_hP9_150_ef_bd-001.{AB}                       hP9       P321                C22       Fe2P                #150      Original Fe2P (part 1)
A3B_hP24_151_3c_2a-001.{AB}                      hP24      P3_{1}12            D0_4      CrCl3               #151      CrCl3 (part 1)      
A2B_hP9_152_c_a-001.{AB}                         hP9       P3_{1}21            -         SiO2                #152      alphaQuartz (part 1)
A3B_hP24_153_3c_2b-001.{AB}                      hP24      P3_{2}12            -         CrCl3               #153      CrCl3 (part 2)      
AB_hP6_154_a_b-001.{AB}                          hP6       P3_{2}21            B9        HgS                 #154      Cinnabar (part 1)   
A3B2_hR5_155_e_c-001.{AB}                        hR5       R32                 D5_e      Ni3S2               #155      Hazelwoodite (part 1)
AB3_hR8_155_c_de-001.{AB}                        hR8       R32                 D0_14     AlF3                #155      AlF3 (part 1)       
AB2_hP9_156_b2c_3a2bc-001.{AB}                   hP9       P3m1                -         CdI2                #156      CdI2 (part 2)       
AB_hP12_156_2ab3c_2ab3c-001.{AB}                 hP12      P3m1                -         CuI                 #156      CuI (part 2)        
AB_hP4_156_ab_ab-001.{AB}                        hP4       P3m1                -         beta-CuI            #156      beta-CuI (part 2)
A3B_hP8_158_d_a-001.{AB}                         hP8       P3c1                -         beta-RuCl3          #158      beta-RuCl3 (part 2) 
A2B3_hP20_159_bc_2c-001.{AB}                     hP20      P31c                -         Bi2O3               #159      Bi2O3 (part 2)      
A4B3_hP28_159_ab2c_2c-001.{AB}                   hP28      P31c                -         alpha-Si3N4         #159      Nierite (part 2)    
A3B_hR4_160_b_a-001.{AB}                         hR4       R3m                 -         H3S                 #160      H3S (part 2)        
A8B5_hR26_160_a3bc_a3b-001.{AB}                  hR26      R3m                 D8_{10}   Cr5Al8              #160      Cr5Al8 (part 2)
AB_hR10_160_5a_5a-001.AB                         hR10      R3m                 B7        SiC                 #160      Moissanite-15R (part 2)
AB_hR6_160_3a_3a-001.AB                          hR6       R3m                 -         CSi                 #160      Moissanite 9R structure (part 1)
AB_hR6_160_b_b-001.{AB}                          hR6       R3m                 B13       NiS                 #160      Millerite (part 1)  
AB2_hP9_162_ad_k-001.{AB}                        hP9       P-31m               -         NV2                 #162      betaV2N (part 1)    
A2B3_hP5_164_d_ad-001.{AB}                       hP5       P-3m1               D5_{2}    La2O3               #164      La2O3 (part 2)      
A3B2_hP5_164_ad_d-001.{AB}                       hP5       P-3m1               D5_13     Al3Ni2              #164      Al3Ni2 (part 1)     
AB2_hP3_164_a_d-001.{AB}                         hP3       P-3m1               C6        CdI2                #164      omega Phase (part 1)
AB2_hP9_164_bd_c2d-001.{AB}                      hP9       P-3m1               -         deltaH^II-NW2       #164      deltaH^II-NW2 (part 2)
A3B_hP24_165_adg_f-001.{AB}                      hP24      P-3c1               -         H3Ho                #165      H3Ho (part 1)       
A3B_hP24_165_bdg_f-001.{AB}                      hP24      P-3c1               D0_{21}   Cu3P                #165      Cu3P (part 2)       
A2B3_hR5_166_c_ac-001.{AB}                       hR5       R-3m                C33       Bi2Te3              #166      Bi2Te3 (part 1)     
A4B3_hR7_166_2c_ac-001.{AB}                      hR7       R-3m                D7_{1}    Al4C3               #166      Al4C3 (part 2)      
A5B2_hR7_166_a2c_c-001.{AB}                      hR7       R-3m                D8_i      B5Mo2               #166      Mo2B5 (part 1)      
A6B_hR7_166_g_a-001.{AB}                         hR7       R-3m                -         CaC6                #166      CaC6 (part 1)       
A7B6_hR13_166_ah_3c-001.{AB}                     hR13      R-3m                D8_5      Fe7W6               #166      Fe7W6 mu-phase (part 1)
AB_hR2_166_a_b-001.AB                            hR2       R-3m                L1_1      CuPt                #166      CuPt (part 1), fcc derivative (Y. Lederer), equivalent to 4 [HTQC] 
A2B3_hR10_167_c_e-001.{AB}                       hR10      R-3c                D5_1      Al2O3               #167      Corundum (part 1)   
A2B3_hR10_167_c_e-002.{AB}                       hR10      R-3c                -         O3V2                #167      O3V2 binary oxide (R. Friedrich), ICSD #94768
A2B3_hP30_169_2a_3a-001.{AB}                     hP30      P6_{1}              -         alpha-Al2S3         #169      alpha-Al2S3 (part 2)
A2B3_hP30_170_2a_3a-001.{AB}                     hP30      P6_{5}              -         Al2S3               #170      Al2S3 (part 2)      
A3B_hP8_173_c_b-001.{AB}                         hP8       P6_{3}              -         PI3                 #173      PI3 (part 2)        
A4B3_hP14_173_bc_c-001.{AB}                      hP14      P6_{3}              -         beta-Si3N4          #173      beta-Si3N4 (part 2) 
A3B2_hP10_176_h_bc-001.{AB}                      hP10      P6_{3}/m            -         Er3Ru2              #176      Er3Ru2 (part 2)
A3B_hP8_176_h_c-001.{AB}                         hP8       P6_{3}/m            -         UCl3                #176      UCl3 (part 2)
A2B_hP36_177_j2lm_n-001.{AB}                     hP36      P622                -         SiO2                #177      SiO2 (part 2)       
AB3_hP24_178_b_ac-001.{AB}                       hP24      P6_{1}22            -         AuF3                #178      AuF3 (part 2)       
AB3_hP24_179_b_ac-001.{AB}                       hP24      P6_{5}22            -         AuF3                #179      AuF3 (part 2)       
A2B_hP18_180_fi_bd-001.{AB}                      hP18      P6_{2}22            C_a       Mg2Ni               #180      Mg2Ni (part 1)      
A2B_hP9_180_j_c-001.{AB}                         hP9       P6_{2}22            C8        SiO2                #180      beta-Quartz (part 1)
AB2_hP9_180_d_j-001.{AB}                         hP9       P6_{2}22            C40       CrSi2               #180      CrSi2 (part 1)      
A2B_hP9_181_j_c-001.{AB}                         hP9       P6_{4}22            -         beta-SiO2           #181      beta-SiO2 (part 2)  
AB3_hP8_182_c_g-001.{AB}                         hP8       P6_{3}22            -         Fe3C                #182      Bainite (part 1)    
AB_hP6_183_c_ab-001.{AB}                         hP6       P6mm                -         CrFe3NiSn5          #183      CrFe3NiSn5 (part 2) 
A3B_hP24_185_ab2c_c-001.{AB}                     hP24      P6_{3}cm            -         Cu3P                #185      Cu3P (part 2)       
A3B_hP8_185_c_a-001.{AB}                         hP8       P6_{3}cm            -         beta-RuCl3          #185      beta-RuCl3 (part 2) 
AB3_hP24_185_c_ab2c-001.{AB}                     hP24      P6_{3}cm            -         Na3As               #185      Na3As (part 2)      
A3B7_hP20_186_c_b2c-001.{AB}                     hP20      P6_{3}mc            D10_{2}   Fe3Th7              #186      Fe3Th7 (part 2)     
AB_hP12_186_a2b_a2b-001.AB                       hP12      P6_{3}mc            B6        SiC                 #186      Moissanite-6H SiC (part 1)
AB_hP4_186_b_a-001.{AB}                          hP4       P6_{3}mc            B12       BN                  #186      Original BN (part 1)
AB_hP4_186_b_b-001.AB                            hP4       P6_{3}mc            B4        ZnS                 #186      Wurtzite (part 1)   
AB_hP8_186_ab_ab-001.AB                          hP8       P6_{3}mc            B5        SiC                 #186      Moissanite-4H SiC (part 1)
AB3_hP4_187_e_fh-001.{AB}                        hP4       P-6m2               -         Re3N                #187      Re3N (part 2)       
AB_hP2_187_d_a-001.AB                            hP2       P-6m2               B_h       WC                  #187      Tungsten Carbide (part 1)
A2B_hP9_189_fg_bc-001.{AB}                       hP9       P-62m               C22       Fe2P                #189      Revised Fe2P (part 1)
A2B_hP18_190_gh_bf-001.{AB}                      hP18      P-62c               -         Li2Sb               #190      Li2Sb (part 2)      
A5B3_hP16_190_bdh_g-001.{AB}                     hP16      P-62c               -         alpha-Sm3Ge5        #190      alpha-Sm3Ge5 (part 2)
AB_hP24_190_i_afh-001.{AB}                       hP24      P-62c               -         FeS                 #190      Troilite (part 2)   
A2B_hP6_191_h_e-001.{AB}                         hP6       P6/mmm              C_h       Cu2Te               #191      Cu2Te (part 1)      
A3B_hP4_191_bc_a-001.{AB}                        hP4       P6/mmm              -         Li3N                #191      Li3 Ni (part 1)     
AB2_hP3_191_a_d-001.{AB}                         hP3       P6/mmm              C32       AlB2                #191      Hexagonal omega (part 1)
AB5_hP6_191_a_cg-001.{AB}                        hP6       P6/mmm              D2_d      CaCu5               #191      CaCu5 (part 1)      
AB_hP6_191_f_ad-001.{AB}                         hP6       P6/mmm              B35       CoSn                #191      CoSn (part 1)       
AB2_hP72_192_m_j2kl-001.{AB}                     hP72      P6/mcc              -         AlPO4               #192      AlPO4 (part 2)      
A5B3_hP16_193_dg_g-001.{AB}                      hP16      P6_{3}/mcm          D8_{8}    Mn5Si3              #193      Mavlyanovite (part 2)
A2B_hP12_194_cg_f-001.{AB}                       hP12      P6_{3}/mmc          C10       SiO2                #194      beta-Tridymite (part 1)
A3B_hP16_194_gh_ac-001.{AB}                      hP16      P6_{3}/mmc          D0_{24}   Ni3Ti               #194      Ni3Ti (part 2)      
A3B_hP8_194_h_c-001.{AB}                         hP8       P6_{3}/mmc          D0_19     Ni3Sn               #194      Ni3Sn (part 1)      
A5B2_hP14_194_abdf_f-001.{AB}                    hP14      P6_{3}/mmc          D8_h      W2B5                #194      W2B5 (part 1)       
A5B2_hP28_194_ahk_ch-001.{AB}                    hP28      P6_{3}/mmc          D8_{11}   Co2Al5              #194      Co2Al5 (part 2)     
AB2_hP12_194_f_ah-001.{AB}                       hP12      P6_{3}/mmc          C14       MgZn2               #194      MgZn2 Hexagonal Laves (part 1)
AB2_hP24_194_ef_fgh-001.{AB}                     hP24      P6_{3}/mmc          C36       MgNi2               #194      MgNi2 Hexagonal Laves (part 1)
AB2_hP6_194_b_f-001.{AB}                         hP6       P6_{3}/mmc          -         CaIn2               #194      CaIn2 (part 1)      
AB2_hP6_194_c_ad-001.{AB}                        hP6       P6_{3}/mmc          B8_2      Ni2In               #194      Ni2In (part 1)      
AB2_hP6_194_c_f-001.{AB}                         hP6       P6_{3}/mmc          C7        MoS2                #194      Molybdenite (part 1)
AB3_hP8_194_c_bf-001.{AB}                        hP8       P6_{3}/mmc          D0_18     Na3As               #194      Na3As (part 1)      
AB_hP12_194_af_bf-001.{AB}                       hP12      P6_{3}/mmc          -         CMo                 #194      CMo structure (part 1)
AB_hP12_194_df_ce-001.{AB}                       hP12      P6_{3}/mmc          B18       CuS                 #194      Covellite (part 1)  
AB_hP4_194_c_a-001.{AB}                          hP4       P6_{3}/mmc          B8_1      NiAs                #194      NiAs (part 1)       
AB_hP4_194_c_d-001.AB                            hP4       P6_{3}/mmc          B_k       BN                  #194      BN (part 1)         
AB_hP8_194_ad_f-001.{AB}                         hP8       P6_{3}/mmc          B_i       AsTi                #194      AsTi (part 1)       
A4B_cI40_197_cde_c-001.{AB}                      cI40      I23                 -         Ga4Ni               #197      Ga4Ni3 (part 1)     
A3B_cP16_198_b_a-001.{AB}                        cP16      P2_{1}3             D1        NH3                 #198      Ammonia (part 1)    
AB_cP8_198_a_a-001.{AB}                          cP8       P2_{1}3             B21       CO                  #198      alphaCO (part 1)    
AB_cP8_198_a_a-002.AB                            cP8       P2_{1}3             B20       FeSi                #198      FeSi (part 1)       
AB_cI16_199_a_a-001.AB                           cI16      I2_{1}3             B_a       CoU                 #199      CoU (part 1)        
A2B11_cP39_200_f_aghij-001.{AB}                  cP39      Pm-3                D8_{c}    Mg2Zn11             #200      Mg2Zn11 (part 2)    
A12B_cI26_204_g_a-001.{AB}                       cI26      Im-3                -         Al12W               #204      Al12W (part 1)      
A3B_cI32_204_g_c-001.{AB}                        cI32      Im-3                D0_2      CoAs3               #204      Skutterudite (part 1)
AB2_cP12_205_a_c-001.{AB}                        cP12      Pa-3                C2        FeS2                #205      Pyrite (part 1)     
AB2_cP12_205_a_c-002.{AB}                        cP12      Pa-3                -         NaO2                #205      binary oxide (R. Friedrich), ICSD #87178, equivalent to 1214 [HTQC]
AB_cP16_205_c_c-001.{AB}                         cP16      Pa-3                -         CuCl                #205      SC16 (part 1)       
A17B15_cP64_207_acfk_eij-001.{AB}                cP64      P432                -         Pd17Se15            #207      Pd17Se15 (part 2)   
A3B_cP16_208_j_b-001.{AB}                        cP16      P4_{2}32            -         PH3                 #208      PH3 (part 2)        
A2B_cI72_211_hi_i-001.{AB}                       cI72      I432                -         SiO2                #211      SiO2 (part 2)       
A2B_cP12_212_c_a-001.{AB}                        cP12      P4_{3}32            -         SrSi2               #212      SrSi2 (part 2)      
A4B9_cP52_215_ei_3efgi-001.{AB}                  cP52      P-43m               D8_{3}    gamma-Cu9Al4        #215      gamma-brass (part 2)
AB4_cP5_215_a_e-001.{AB}                         cP5       P-43m               -         Fe4C                #215      Fe4C (part 1)       
AB5_cF24_216_a_ce-001.{AB}                       cF24      F-43m               C15_b     AuBe5               #216      AuBe5 (part 1)      
AB_cF8_216_c_a.AB                                cF8       F-43m               B3        ZnS                 #216      ZnS (part 1)       
A4B_cI10_217_c_a-001.{AB}                        cI10      I-43m               -         SiF4                #217      SiF4 (part 1)       
A5B8_cI52_217_ce_cg-001.{AB}                     cI52      I-43m               -         Cu5Zn8              #217      gammaBrass (part 1) 
A15B4_cI76_220_ae_c-001.{AB}                     cI76      I-43d               D8_{6}    Cu15Si4             #220      Cu15Si4 (part 2)    
A3B2_cI40_220_d_c-001.{AB}                       cI40      I-43d               D5_c      Pu2C3               #220      Pu2C3 (part 1)      
A4B3_cI28_220_c_a-001.{AB}                       cI28      I-43d               D7_{3}    Th3P4               #220      Th3P4 (part 2)      
A3B_cP4_221_d_a.{AB}                             cP4       Pm-3m               D0_9      ReO3                #221      alpha-ReO3 (part 1)  
A6B_cP7_221_f_a-001.{AB}                         cP7       Pm-3m               D2_1      CaB6                #221      CaB6 (part 1)       
AB11_cP36_221_c_agij-001.{AB}                    cP36      Pm-3m               D2_e      BaHg11              #221      BaHg11 (part 1)     
AB3_cP4_221_a_c.{AB}                             cP4       Pm-3m               L1_2      Cu3Au               #221      Cu3Au (part 1), fcc derivative (Y. Lederer), equivalent to 25,26 [HTQC] 
AB_cP2_221_b_a.AB                                cP2       Pm-3m               B2/G0_{8} CsCl/(NH4)(NO3)     #221      CsCl (part 1), equivalent to 61,1026.AB,1205.AB [HTQC], equivalent to NH4NO3 I (G0_{8}, http://aflow.org/prototype-encyclopedia/AB_cP2_221_a_b.NH4.NO3.html, part 3)
AB_cP6_221_c_d.AB                                cP6       Pm-3m               -         NbO                 #221      NbO (part 1)        
A3B_cP8_223_c_a.{AB}                             cP8       Pm-3n               A15       Cr3Si               #223      Cr3Si (part 1)      
A2B_cP6_224_b_a.{AB}                             cP6       Pn-3m               C3        Cu2O                #224      Cuprite (part 1)    
A12B_cF52_225_i_a-001.{AB}                       cF52      Fm-3m               D2_f      UB12                #225      UB12 (part 1)       
A23B6_cF116_225_bd2f_e-001.{AB}                  cF116     Fm-3m               D8_{a}    Th6Mn23             #225      Th6Mn23 (part 2)    
A6B23_cF116_225_e_acfh-001.{AB}                  cF116     Fm-3m               D8_4      Cr23C6              #225      Cr23C6 (part 1)     
A7B_cF32_225_bd_a.{AB}                           cF32      Fm-3m               -/L1_{a}  Ca7Ge/CuPt3         #225      Ca7Ge/CuPt3 (part 1/part 3)      
AB2_cF12_225_a_c.{AB}                            cF12      Fm-3m               C1        CaF2                #225      CaF2 (part 1)
AB3_cF16_225_a_bc.{AB}                           cF16      Fm-3m               D0_3      BiF3                #225      BiF3 (part 1)       
AB_cF8_225_a_b.AB                                cF8       Fm-3m               B1        NaCl                #225      Rock Salt (part 1), equivalent to 201,720.AB,720.BA,1009.AB,1200.AB [HTQC] 
AB13_cF112_226_a_bi-001.{AB}                     cF112     Fm-3c               D2_{3}    NaZn13              #226      NaZn13 (part 2)     
A2B_cF24_227_c_a.{AB}                            cF24      Fd-3m               C9        SiO2                #227      Ideal beta-Cristobalite (part 1)
A2B_cF24_227_d_a.{AB}                            cF24      Fd-3m               C15       Cu2Mg               #227      Cu2Mg Cubic Laves (part 1)
A3B4_cF56_227_ad_e-001.{AB}                      cF56      Fd-3m               D7_{2}    Co3O4               #227      Spinel (part 2)     
AB2_cF48_227_c_e-001.{AB}                        cF48      Fd-3m               -         CTi2                #227      CTi2 (part 1)       
AB2_cF96_227_e_cf-001.{AB}                       cF96      Fd-3m               -         NiTi2               #227      NiTi2 (part 1)      
AB_cF16_227_a_b.AB                               cF16      Fd-3m               B32       NaTl                #227      NaTl (part 1)       
A6B_cF224_228_h_c-001.{AB}                       cF224     Fd-3c               -         TeO6H6              #228      TeO6H6 (part 2)     
A2B7_cI54_229_e_afh-001.{AB}                     cI54      Im-3m               L2_2      Sb2Tl7              #229      L22 (part 1)        
A3B10_cI52_229_e_fh-001.{AB}                     cI52      Im-3m               D8_{1}    gamma-Fe3Zn10       #229      gamma-brass (part 2)
A3B_cI8_229_b_a.{AB}                             cI8       Im-3m               -         H3S                 #229      High Pressure H3S (part 1)
A4B3_cI14_229_c_b.{AB}                           cI14      Im-3m               -         Pt3O4               #229      Pt3O4 (part 1)      
A4B_cI10_229_c_a.{AB}                            cI10      Im-3m               -         beta-Hg4Pt          #229      beta-Hg4Pt (part 2) 
A7B3_cI40_229_df_e-001.{AB}                      cI40      Im-3m               D8_{f}    Ir3Ge7              #229      Ir3Ge7 (part 2)     
A4B3_cI112_230_af_g-001.{AB}                     cI112     Ia-3d               -         Ga4Ni3              #230      Ga4Ni3 (part 1)     

*****************************************************************************************************************************************************
ANRL - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC2_aP16_1_4a_4a_8a-001.{ABC}                   aP16      P1                  -         AsKSe2              #1        AsKSe2 (part 1)     
A3BC_mC10_8_ab_a_a-001.{ABC}                     mC10      Cm                  -         Pb(Zr0.52Ti0.48)O3  #8        Monoclinic PZT [PbO3] (part 1)
ABC2_mP8_10_ac_eh_mn-001.ABC/BCA/CAB             mP8       P2/m                -         AuAgTe2             #10       Muthmannite (part 2)
ABC3_mP10_11_e_e_ef-001.{ABC}                    mP10      P2_{1}/m            G0_6      KClO3               #11       KClO3 (part 1)      
A7B2C2_mC22_12_aij_h_i-001.{ABC}                 mC22      C2/m                S2_{1}    [Sc,Y]2Si2O7        #12       Thortveitite (part 2)
ABC4_mP12_13_e_a_2g-001.{ABC}                    mP12      P2/c                E1_b      AgAuTe4             #13       Sylvanite (part 1)  
ABC4_oP12_16_ag_cd_2u-001.{ABC}                  oP12      P222                -         AlPS4               #16       AlPS4 (part 1)      
A2BC2_oF40_22_fi_ad_gh-001.{ABC}                 oF40      F222                -         CeRu2B2             #22       CeRu2B2 (part 2)    
ABC2_oI16_23_ab_i_k-001.{ABC}                    oI16      I222                -         NaFeS2              #23       NaFeS2 (part 2)     
ABC4_oI12_23_a_b_k-001.{ABC}                     oI12      I222                -         BPS4                #23       BPS4 (part 2)       
ABC_oP12_29_a_a_a-001.ABC/BCA/CAB                oP12      Pca2_{1}            -         CoAsS               #29       Cobaltite (part 2)  
A5B3C15_oP46_30_a2c_bc_a7c-001.{ABC}             oP46      Pnc2                -         Bi5Nb3O15           #30       Bi5Nb3O15 (part 2)  
ABC3_oP20_30_2a_c_3c-001.{ABC}                   oP20      Pnc2                -         CuBrSe3             #30       CuBrSe3 (part 2)    
AB3C4_oP16_31_a_ab_2ab-001.{ABC}                 oP16      Pmn2_{1}            H2_5      AsCu3S4             #31       Enargite (part 1)   
A13B2C2_oP34_32_a6c_c_c-001.{ABC}                oP34      Pba2                -         Re2O5[SO4]2         #32       Re2O52 (part 2)     
AB3C4_oP32_33_a_3a_4a-001.{ABC}                  oP32      Pna2_{1}            -         AsK3S4              #33       AsK3S4 (part 1)     
A2B8C_oP22_34_c_4c_a-001.{ABC}                   oP22      Pnn2                -         TiAl2Br8            #34       TiAl2Br8 (part 2)   
AB8C2_oC22_35_a_ab3e_e-001.{ABC}                 oC22      Cmm2                -         V2MoO8              #35       V2MoO8 (part 2)     
A2B5C2_oC36_37_d_c2d_d-001.{ABC}                 oC36      Ccc2                -         Li2Si2O5            #37       Li2Si2O5 (part 2)   
A2BC_oC8_38_e_a_b-001.{ABC}                      oC8       Amm2                -         C2CeNi              #38       C2CeNi (part 1)     
A9BC_oC44_39_3c3d_a_c-001.{ABC}                  oC44      Aem2                -         VPCl9               #39       VPCl9 (part 2)      
AB2C_oC16_40_a_2b_b-001.{ABC}                    oC16      Ama2                -         K2CdPb              #40       K2CdPb (part 2)     
A2BC2_oI20_45_c_b_c-001.ABC/BCA/CAB              oI20      Iba2                -         MnGa2Sb2            #45       MnGa2Sb2 (part 2)   
ABC_oI36_46_ac_bc_3b-001.{ABC}                   oI36      Ima2                -         TiFeSi              #46       TiFeSi (part 2)     
A2BC4_oP28_50_ij_ac_ijm-001.{ABC}                oP28      Pban                -         La2NiO4             #50       La2NiO4 (part 2)    
A3BC2_oP48_50_3m_m_2m-001.{ABC}                  oP48      Pban                -         alpha-Tl2TeO3       #50       alpha-Tl2TeO3 (part 2)
ABC2_oP16_53_h_e_gh-001.{ABC}                    oP16      Pmna                -         TaNiTe2             #53       TaNiTe2 (part 2)    
ABC3_oP20_53_e_g_hi-001.{ABC}                    oP20      Pmna                -         CuBrSe3             #53       CuBrSe3 (part 2)    
ABC3_oP20_54_e_d_cf-001.{ABC}                    oP20      Pcca                -         BiGaO3              #54       BiGaO3 (part 2)     
ABC_oP6_59_a_a_a-001.{ABC}                       oP6       Pmmn                -         CNCl                #59       CNCl (part 1)       
ABC_oP6_59_a_b_a-001.{ABC}                       oP6       Pmmn                E0_{5}    FeOCl               #59       FeOCl (part 2)
A2B4C_oP28_62_ac_2cd_c-001.{ABC}                 oP28      Pnma                S1_{2}    Mg2SiO4             #62       Forsterite (part 2) 
A2B4C_oP28_62_ac_2cd_c-002.{ABC}                 oP28      Pnma                -         Ag2CrO4             #62       Ag2CrO4 ternary oxide (R. Friedrich), ICSD #16298
AB2C3_oP24_62_c_d_cd-001.{ABC}                   oP24      Pnma                E9_{e}    CuFe2S3             #62       Cubanite (part 2)   
AB2C_oP16_62_c_2c_c-001.{ABC}                    oP16      Pnma                F5_6      CuSbS2              #62       Chalcostibite (part 1)
AB3C_oP20_62_c_cd_a-001.{ABC}                    oP20      Pnma                -         CaTiO3              #62       CaTiO3 Pnma Perovskite (part 1)
AB3C_oP20_62_c_cd_a-002.{ABC}                    oP20      Pnma                -         CaO3Zr              #62       CaO3Zr ternary oxide (R. Friedrich), ICSD #97463
AB4C_oP24_62_c_2cd_c-001.{ABC}                   oP24      Pnma                H0_{2}    BaSO4               #62       Barite (part 2)     
A2BC3_oC24_63_e_c_cg-001.{ABC}                   oC24      Cmcm                -         KFe2S3              #63       Rasvumite (part 2)  
A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h-001.{ABC}    oC260     Cmcm                -         La43Ni17Mg5         #63       La43Ni17Mg5 (part 2)
AB2C_oC16_63_c_2c_c-001.{ABC}                    oC16      Cmcm                -         SrCuO2              #63       SrCuO2 (part 1)     
AB3C_oC20_63_a_cf_c-001.{ABC}                    oC20      Cmcm                -         MgSiO3/AsCV3        #63       Post-perovskite/V3AsC (part 2/part 3); V3AsC is structurally similar to the post-perovskite
AB4C_oC24_63_a_fg_c-001.{ABC}                    oC24      Cmcm                -         MgO4S               #63       MgSO4 (part 2)      
AB4C_oC24_63_c_fg_c-001.{ABC}                    oC24      Cmcm                H0_{1}    CaSO4               #63       Anhydrite (part 2)  
A2B2C_oC80_64_efg_efg_df-001.{ABC}               oC80      Cmce                -         MgB2C2              #64       MgB2C2 (part 1)     
A2B4C_oC28_66_l_kl_a-001.{ABC}                   oC28      Cccm                -         SrAl2Se4            #66       SrAl2Se4 (part 2)   
A2BC_oC16_67_ag_b_g-001.{ABC}                    oC16      Cmme                -         Al2CuIr             #67       Al2CuIr (part 2)    
ABC2_oC16_67_b_g_ag-001.{ABC}                    oC16      Cmme                -         HoCuP2              #67       HoCuP2 (part 2)     
ABC_oI12_71_h_j_g-001.{ABC}                      oI12      Immm                -         NbPS                #71       NbPS (part 2)       
A2BC_tP16_76_2a_a_a-001.{ABC}                    tP16      P4_{1}              -         LaRhC2              #76       LaRhC2 (part 2)     
A2B7C2_tP88_78_4a_14a_4a-001.{ABC}               tP88      P4_{3}              -         Sr2As2O7            #78       Sr2As2O7 (part 2)   
A2BC2_tI20_79_c_2a_c-001.{ABC}                   tI20      I4                  -         TlZn2Sb2            #79       TlZn2Sb2 (part 2)   
A2BC4_tI14_82_bc_a_g-001.{ABC}                   tI14      I-4                 E3        CdAl2S4             #82       CdAl2S4 (part 1)    
AB4C_tI12_82_c_g_a-001.ABC/BCA/CAB               tI12      I-4                 H0_7      BPO4                #82       BPO4 (part 1)       
ABC_tP24_91_d_d_d-001.{ABC}                      tP24      P4_{1}22            -         ThBC                #91       ThBC (part 2)       
A14B3C5_tP44_94_c3g_ad_bg-001.{ABC}              tP44      P4_{2}2_{1}2        -         Na5Fe3F14           #94       Na5Fe3F14 (part 2)  
A6B2C_tP18_94_eg_c_a-001.{ABC}                   tP18      P4_{2}2_{1}2        -         Li2MoF6             #94       Li2MoF6 (part 2)    
ABC_tP24_95_d_d_d-001.{ABC}                      tP24      P4_{3}22            -         ThBC                #95       ThBC (part 2)       
AB8C2_tI44_97_e_2k_cd-001.{ABC}                  tI44      I422                -         Ta2Se8I             #97       Ta2Se8I (part 2)    
A3BC_tP5_99_bc_a_b-001.{ABC}                     tP5       P4mm                -         Pb(Zr0.52Ti0.48)O3  #99       Tetragonal PZT [PbO3] (part 1)
A3B11C6_tP40_100_ac_bc2d_cd-001.{ABC}            tP40      P4bm                -         Ce3Si6N11           #100      Ce3Si6N11 (part 2)  
A7B7C2_tP32_101_bde_ade_d-001.{ABC}              tP32      P4_{2}cm            -         gamma-MgNiSn        #101      gamma-MgNiSn (part 2)
A5B5C4_tP28_104_ac_ac_c-001.{ABC}                tP28      P4nc                -         Ba5In4Bi5           #104      Ba5In4Bi5 (part 2)  
AB6C4_tP22_104_a_2ac_c-001.{ABC}                 tP22      P4nc                -         Tl4HgI6             #104      Tl4HgI6 (part 2)    
A2BC2_tP20_105_f_ac_2e-001.{ABC}                 tP20      P4_{2}mc            -         BaGe2As2            #105      BaGe2As2 (part 2)   
ABC_tI12_109_a_a_a-001.ABC/BCA/CAB               tI12      I4_{1}md            -         LaPtSi              #109      LaPtSi (part 2)     
A2BC8_tI176_110_2b_b_8b-001.{ABC}                tI176     I4_{1}cd            -         Be[BH4]2            #110      Be2 (part 2)        
AB4C_tP12_112_b_n_e-001.{ABC}                    tP12      P-42c               -         alpha-CuAlCl4       #112      alpha-CuAlCl4 (part 2)
A4BC2_tI28_120_i_d_e-001.{ABC}                   tI28      I-4c2               -         KAu4Sn2             #120      KAu4Sn2 (part 2)    
ABC2_tI16_122_a_b_d-001.{ABC}                    tI16      I-42d               E1_1      CuFeS2              #122      Chalcopyrite (part 1)
AB5C_tP7_123_b_ci_a-001.{ABC}                    tP7       P4/mmm              -         HoCoGa5             #123      HoCoGa5 (part 1)    
ABC2_tP4_123_d_a_f-001.{ABC}                     tP4       P4/mmm              -         CaCuO2              #123      CaCuO2 (part 1)     
AB4C_tP12_124_a_m_c-001.ABC/BCA/CAB              tP12      P4/mcc              -         Nb4CoSi             #124      Nb4CoSi (part 2)    
ABC4_tP12_125_a_b_m-001.{ABC}                    tP12      P4/nbm              -         KCeSe4              #125      KCeSe4 (part 2)     
A2BC4_tP28_126_cd_e_k-001.{ABC}                  tP28      P4/nnc              -         BiAl2S4             #126      BiAl2S4 (part 2)    
A6B2C_tP18_128_eh_d_a-001.{ABC}                  tP18      P4/mnc              -         K2SnCl6             #128      K2SnCl6 (part 2)    
A7B2C_tP40_128_egi_h_e-001.{ABC}                 tP40      P4/mnc              E9_{a}    FeCu2Al7            #128      FeCu2Al7 (part 2)   
ABC_tP6_129_c_a_c-001.{ABC}                      tP6       P4/nmm              E0_1      PbFCl               #129      Matlockite (part 1) 
A2BC4_tP28_130_f_c_g-001.{ABC}                   tP28      P4/ncc              -         CuBi2O4             #130      CuBi2O4 (part 2)    
AB6C_tP16_132_d_io_a-001.{ABC}                   tP16      P4_{2}/mcm          -         AgUF6               #132      AgUF6 (part 2)      
A4B2C_tP28_135_gh_h_d-001.{ABC}                  tP28      P4_{2}/mbc          -         ZnSb2O4             #135      ZnSb2O4 (part 2)    
A4BC4_tP18_137_g_b_g-001.{ABC}                   tP18      P4_{2}/nmc          -         CeCo4B4             #137      CeCo4B4 (part 2)    
AB2C4_tI14_139_a_e_ce-001.{ABC}                  tI14      I4/mmm              -         (La,Ba)2CuO4        #139      0201 [(La,Ba)2CuO4] (part 1)
AB2C4_tI14_139_a_e_ce-002.{ABC}                  tI14      I4/mmm              -         O4Sr2Ti             #139      O4Sr2Ti ternary oxide (R. Friedrich), ICSD #157402, equivalent to K2NiF4 (http://aflow.org/prototype-encyclopedia/A4B2C_tI14_139_ce_e_a.html, part 3)
A4BC_tI24_141_h_b_a-001.{ABC}                    tI24      I4_{1}/amd          -         ZrSiO4              #141      Zircon (part 1)     
ABC4_tI96_142_e_ab_2g-001.{ABC}                  tI96      I4_{1}/acd          -         PPrS4               #142      PPrS4 (part 1)      
A4B14C3_hP21_143_bd_ac4d_d-001.{ABC}             hP21      P3                  -         Ta3Al4O13[OH]       #143      Simpsonite (part 2) 
A4B6C_hP11_143_bd_2d_a-001.{ABC}                 hP11      P3                  -         ScRh6P4             #143      ScRh6P4 (part 2)    
A3BC_hR5_146_b_a_a-001.{ABC}                     hR5       R3                  -         gamma-Ag3SI         #146      gamma-Ag3SI (part 2)
ABC3_hR10_146_2a_2a_2b-001.{ABC}                 hR10      R3                  -         FePSe3              #146      FePSe3 (part 2)     
A2B4C_hR42_148_2f_4f_f-001.{ABC}                 hR42      R-3                 S1_{3}    Be2SiO4             #148      Phenakite (part 2)  
AB3C_hR10_148_c_f_c-001.ABC/BCA/CAB              hR10      R-3                 -         FeTiO3              #148      Ilmenite (part 1)   
A5B6C2_hP13_157_2ac_2c_b-001.{ABC}               hP13      P31m                -         Ag5Pb2O6            #157      Ag5Pb2O6 (part 2)   
ABC_hR3_160_a_a_a-001.ABC/BAC/CAB                hR3       R3m                 F0_{2}    COS                 #160      Carbonyl Sulphide (part 2)
ABC3_hR10_161_a_a_b-001.{ABC}                    hR10      R3c                 -         LiNbO3              #161      Ferroelectric LiNbO3 (part 1)
ABC2_hP4_164_a_b_d-001.{ABC/BCA/CAB              hP4       P-3m1               -         CuNiSb2             #164      CuNiSb2 (part 2)    
ABC2_hR4_166_a_b_c-001.{ABC}                     hR4       R-3m                F5_1      CrNiS2              #166      Caswellsilverite (part 1)
ABC2_hR4_166_a_b_c-002.ABC/BCA/CAB               hR4       R-3m                -         ABC2                #166      rocksalt derivative (Y. Lederer)
ABC2_hR4_166_a_b_c-003.ABC/BCA/CAB               hR4       R-3m                -         AlLiO2              #166      AlLiO2 ternary oxide (R. Friedrich), ICSD #28288 
ABC2_hR4_166_a_b_c-004.{ABC}                     hR4       R-3m                -         CuFeO2              #166      CuFeO2 ternary oxide (R. Friedrich), ICSD #246912, equivalent to rhombohedral delafossite (CuFeO2, http://aflow.org/prototype-encyclopedia/ABC2_hR4_166_a_b_c.CuFeO2.html, part 3)
ABC_hR6_166_c_c_c-001.{ABC}                      hR6       R-3m                -         SmSI                #166      SmSI (part 2)       
AB3C_hR10_167_b_e_a-001.{ABC}                    hR10      R-3c                -         PrNiO3              #167      PrNiO3 (part 2)     
ABC2_hR24_167_e_e_2e-001.{ABC}                   hR24      R-3c                F5_{13}   KBO2                #167      KBO2 (part 2)       
ABC3_hR10_167_a_b_e-001.{ABC}                    hR10      R-3c                -         LiNbO3              #167      Paraelectric LiNbO3 (part 1)
ABC3_hR10_167_a_b_e-002.{ABC}                    hR10      R-3c                G0_1      CaCO3               #167      Calcite (part 1)    
A2B13C4_hP57_168_d_c6d_2d-001.{ABC}              hP57      P6                  -         K2Ta4O9F4           #168      K2Ta4O9F4 (part 2)  
AB4C_hP72_168_2d_8d_2d-001.{ABC}                 hP72      P6                  -         Al[PO4]             #168      Al (part 2)         
A10B2C_hP39_171_5c_c_a-001.{ABC}                 hP39      P6_{2}              -         Sr[S2O6][H2O]4      #171      Sr4 (part 2)        
A10B2C_hP39_172_5c_c_a-001.{ABC}                 hP39      P6_{4}              -         Sr[S2O6][H2O]4      #172      Sr4 (part 2)        
A12B7C2_hP21_174_2j2k_ajk_cf-001.{ABC}           hP21      P-6                 -         Fe12Zr2P7           #174      Fe12Zr2P7 (part 2)  
ABC_hP12_174_cj_fk_aj-001.ABC/BCA/CAB            hP12      P-6                 -         GdSI                #174      GdSI (part 2)       
A8B7C6_hP21_175_ck_aj_k-001.{ABC}                hP21      P6/m                -         Nb7Ru6B8            #175      Nb7Ru6B8 (part 2)   
ABC_hP36_175_jk_jk_jk-001.{ABC}                  hP36      P6/m                -         Mg[NH]              #175      Mg (part 2)         
A3B3C_hP14_176_h_h_c-001.{ABC}                   hP14      P6_{3}/m            -         Fe3Te3Tl            #176      Fe3Te3Tl (part 3)
ABC_hP3_183_a_a_a-001.ABC/BCA/CAB                hP3       P6mm                -         AuCN                #183      AuCN (part 2)       
AB4C_hP72_184_d_4d_d-001.ABC/BCA/CAB             hP72      P6cc                -         Al[PO4]             #184      Al (part 2)         
A3BC_hP30_185_cd_c_ab-001.{ABC}                  hP30      P6_{3}cm            -         KNiCl3              #185      KNiCl3 (part 2)     
A5B3C_hP18_186_2a3b_2ab_b-001.{ABC}              hP18      P6_{3}mc            E9_4      Al5C3N              #186      Al5C3N (part 1)     
ABC_hP3_187_a_d_f-001.ABC/BCA/CAB                hP3       P-6m2               -         BaPtSb              #187      BaPtSb (part 1)     
A3BC_hP10_188_k_c_a-001.ABC/BCA/CAB              hP10      P-6c2               -         LiScI3              #188      LiScI3 (part 2)
AB9C4_hP28_188_e_kl_ak-001.{ABC}                 hP28      P-6c2               S3_{2}    BaSi4O9             #188      BaSi4O9 (part 2)
AB4C_hP6_191_a_h_b-001.ABC/BCA/CAB               hP6       P6/mmm              -         AlB4Mg              #191      AlB4Mg (part 1)     
A9B3C_hP26_194_hk_h_a-001.{ABC}                  hP26      P6_{3}/mmc          E9_{c}    Al9Mn3Si            #194      Al9Mn3Si (part 2)   
AB3C4_hP16_194_c_af_ef-001.{ABC}                 hP16      P6_{3}/mmc          -         AlN3Ti4             #194      AlN3Ti4 (part 1)    
ABC2_hP8_194_d_a_f-001.{ABC}                     hP8       P6_{3}/mmc          -         AlCCr2              #194      AlCCr2 (part 1)     
ABC_hP6_194_c_d_a-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         LiBC                #194      LiBC (part 1)       
A12BC4_cP34_195_2j_ab_2e-001.{ABC}               cP34      P23                 -         PrRu4P12            #195      PrRu4P12 (part 2)   
A12B2C_cF60_196_h_bc_a-001.{ABC}                 cF60      F23                 -         Cu2Fe[CN]6          #196      Cu2Fe6 (part 2)     
ABC3_cP20_198_a_a_b-001.{ABC}                    cP20      P2_{1}3             G3        NaClO3              #198      Sodium Chlorate (part 2)
ABC_cP12_198_a_a_a-001.{ABC}                     cP12      P2_{1}3             F0_1      NiSSb               #198      Ullmanite (part 1)  
AB3C_cP60_201_be_fh_g-001.{ABC}                  cP60      Pn-3                -         KSbO3               #201      KSbO3 (part 2)
A6B6C_cF104_202_h_h_c-001.{ABC}                  cF104     Fm-3                -         KB6H6               #202      KB6H6 (part 2)      
AB3C16_cF160_203_a_bc_eg-001.{ABC}               cF160     Fd-3                -         Rb3AsSe16           #203      Rb3AsSe16 (part 2)
AB32C48_cI162_204_a_2efg_2gh-001.{ABC}           cI162     Im-3                -         Mg32(Al,Zn)49       #204      Bergman [Mg32(Al,Zn)49] (part 1)
A2B3C6_cP264_205_2d_ab2c2d_6d-001.{ABC}          cP264     Pa-3                -         Ca3Al2O6            #205      Ca3Al2O6 (part 2)   
AB3C2_cI96_206_c_e_ad-001.{ABC}                  cI96      Ia-3                E9_{d}    AlLi3N2             #206      AlLi3N2 (part 2)    
AB3C6_cI80_206_b_d_e-001.{ABC}                   cI80      Ia-3                D5_3      (Mn,Fe)2O3          #206      Bixbyite (part 1)
A24BC_cF104_209_j_a_b-001.{ABC}                  cF104     F432                -         F6KP                #209      F6KP (part 2)       
A12B6C_cF608_210_4h_2h_e-001.{ABC}               cF608     F4_{1}32            -         Te[OH]6             #210      Te6 (part 2)        
A3B3C_cI56_214_g_h_a-001.{ABC}                   cI56      I4_{1}32            -         Ca3PI3              #214      Ca3PI3 (part 2)     
A3BC2_cI48_214_f_a_e-001.{ABC}                   cI48      I4_{1}32            -         Ag3AuTe2            #214      Petzite (part 2)    
A3B4C_cP8_215_d_e_a-001.{ABC}                    cP8       P-43m               H2_4      Cu3S4V              #215      Sulvanite (part 1)  
AB3C4_cP8_215_a_c_e-001.{ABC}                    cP8       P-43m               -         Cu3AsS4             #215      Cubic Lazarevicite (part 1)
ABC_cF12_216_b_c_a.ABC/BCA/CAB                   cF12      F-43m               C1_b      AgAsMg              #216      Half-Heusler (part 1)
A3B4C_cP16_218_c_e_a-001.{ABC}                   cP16      P-43n               H2_{1}    Ag3[PO4]            #218      Ag3 (part 2)
A2B3C6_cP33_221_cd_ag_fh-001.{ABC}               cP33      Pm-3m               E9_{1}    Ca3Al2O6            #221      Ca3Al2O6 (part 2)   
AB3C_cP5_221_a_c_b.{ABC}                         cP5       Pm-3m               E2_1      CaTiO3              #221      Cubic Perovskite (part 1)
A5B3C16_cP96_222_ce_d_fi-001.{ABC}               cP96      Pn-3n               -         Ce5Mo3O16           #222      Ce5Mo3O16 (part 2)  
A6B2C_cF36_225_e_c_a-001.{ABC}                   cF36      Fm-3m               J1_{1}    K2PtCl6             #225      K2PtCl6 (part 2)    
A9B16C7_cF128_225_acd_2f_be-001.{ABC}            cF128     Fm-3m               -         Cr9Fe16Ni7          #225      Model of Ferrite (part 1)
AB18C8_cF108_225_a_eh_f-001.{ABC}                cF108     Fm-3m               -         CrFe18Ni8           #225      Model of Austenite (part 1)
AB2C_cF16_225_a_c_b.ABC/BCA/CAB                  cF16      Fm-3m               L2_1      AlCu2Mn             #225      Heusler (part 1), B2 derivative (Y. Lederer), equivalent to T0001 [HTQC] 
A2B2C7_cF88_227_c_d_af-001.{ABC}                 cF88      Fd-3m               -         Eu2Ir2O7            #227      Pyrochlore Iridate (part 2)
A2BC4_cF56_227_d_a_e-001.{ABC}                   cF56      Fd-3m               H1_1      Al2MgO4             #227      Spinel (part 1)     
AB3C3_cF112_227_c_de_f-001.{ABC}                 cF112     Fd-3m               E9_3      Fe3W3C              #227      Fe3W3C (part 1)     
AB12C3_cI32_229_a_h_b-001.{ABC}                  cI32      Im-3m               -         CrFe12Ni3           #229      Model of Austenite (part 1)
AB4C3_cI16_229_a_c_b.{ABC}                       cI16      Im-3m               -         CrFe4Ni3            #229      Model of Ferrite (part 1)

*****************************************************************************************************************************************************
ANRL - QUATERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC6D2_mC40_15_e_e_3f_f-001.{ABCD}               mC40      C2/c                -/S4_{1}  CaFeO6Si2/CaMgO6Si2 #15       Esseneite/Diopside (part 1/part 3), Diopside (S4_{1}) is structurally similar to esseneite 
AB7CD2_oI44_24_a_b3d_c_ac-001.{ABCD}             oI44      I2_{1}2_{1}2_{1}    -         Na2MgAlF7           #24       Weberite (part 2)   
A6B4C16D_oP108_27_abcd4e_4e_16e_e-001.{ABCD}     oP108     Pcc2                -         Ca4Al6O16S          #27       Ca4Al6O16S (part 2) 
A2B3C7D_oP13_47_t_aq_eqrs_h-001.{ABCD}           oP13      Pmmm                -         YBa2Cu3O7-x         #47       1212C [YBa2Cu3O7-x] (part 1)
A2B8CD_oP24_48_k_2m_d_b-001.{ABCD}               oP24      Pnnn                -         alpha-RbPr[MoO4]2   #48       alpha-RbPr2 (part 2)
AB2C8D_oP24_49_g_q_2qr_e-001.{ABCD}              oP24      Pccm                -         CsPr[MoO4]2         #49       CsPr2 (part 2)      
ABCD_oP16_57_d_c_d_d-001.{ABCD}                  oP16      Pbcm                F5_9      KCNS                #57       KCNS (part 1)       
ABCD3_oI48_73_d_e_e_ef-001.{ABCD}                oI48      Ibca                -         KAg[CO3]            #73       KAg (part 2)        
AB2C12D4_tP76_75_2a2b_2d_12d_4d-001.{ABCD}       tP76      P4                  -         BaCr2Ru4O12         #75       BaCr2Ru4O12 (part 2)
A2B6CD7_tP64_77_2d_6d_d_ab6d-001.{ABCD}          tP64      P4_{2}              -         MgB2O(OH)6          #77       Pinnoite (part 2)   
A17BC4D_tP184_89_17p_p_4p_io-001.{ABCD}          tP184     P422                -         C17FeO4Pt           #89       C17FeO4Pt (part 2)  
A4B2C13D_tP40_90_g_d_cef2g_c-001.{ABCD}          tP40      P42_{1}2            -         Na4Ti2Si8O22[H2O]4  #90       Na4Ti2Si8O224 (part 2)
A2B8CD_tI24_97_d_k_a_b-001.ABCD/ACDB/ADBC/       tI24      I422                -         NaGdCu2F8           #97       NaGdCu2F8 (part 2)  
BADC/BCAD/BDCA/CABD/CBDA/CDAB/DACB/DBAC/DCBA
A2B8C2D_tP26_100_c_abcd_c_a-001.{ABCD}           tP26      P4bm                -         Ba2TiSi2O8          #100      Fresnoite (part 2)  
A3BC3D_tP64_106_3c_c_3c_c-001.{ABCD}             tP64      P4_{2}bc            -         NaZn[OH]3           #106      NaZn3 (part 2)      
A2BC7D2_tP24_113_e_a_cef_e-001.{ABCD}            tP24      P-42_{1}m           S5_{3}    Ca2MgSi2O7          #113      Akermanite (part 2) 
A2BC4D_tI16_121_d_a_i_b-001.ABCD/ACDB/ADBC/      tI16      I-42m               H2_6      Cu2FeS4Sn           #121      Stannite (part 1)   
BADC/BCAD/BDCA/CABD/CBDA/CDAB/DACB/DBAC/DCBA
A4BC4D_tP10_123_gh_a_i_d-001.{ABCD}              tP10      P4/mmm              -         CaRbFe4As4          #123      CaRbFe4As4 (part 2)
ABCD_tP8_129_c_b_a_c-001.{ABCD}                  tP8       P4/nmm              -         AsCuSiZr            #129      AsCuSiZr (part 1)   
A2B2C4D_tP18_132_e_i_o_d-001.{ABCD}              tP18      P4_{2}/mcm          -         Rb2TiCu2Se4         #132      Rb2TiCu2S4 (part 2) 
AB4C7D_hP26_159_b_ac_a2c_b-001.{ABCD}            hP26      P31c                -         YbBaCo4O7           #159      YbBaCo4O7 (part 2)  
AB2CD2_hP36_163_h_i_bf_i-001.{ABCD}              hP36      P-31c               F5_10     KAg(CN)2            #163      KAg2 (part 1)       
A8BC3D6_hP18_189_bfh_a_g_i-001.{ABCD}            hP18      P-62m               E9_{b}    pi-FeMg3Al8Si6      #189      pi-FeMg3Al8Si6 (part 2)
A9BC3D5_hP18_189_fi_a_g_bh-001.{ABCD}            hP18      P-62m               -         pi-FeMg3Al9Si5      #189      pi-FeMg3Al9Si5 (part 2)
A2B3C18D6_hP58_192_c_f_lm_l-001.{ABCD}           hP58      P6/mcc              G3_{1}    Be3Al2Si6O18        #192      Beryl (part 2)      
A12B36CD12_cF488_210_h_3h_a_fg-001.{ABCD}        cF488     F4_{1}32            -         MgB12H12[H2O]12     #210      MgB12H1212 (part 2)
ABCD_cF16_216_c_d_b_a.ABCD/BCAD/CABD             cF16      F-43m               -         LiMgAuSn            #216      Quaternary Heusler (part 2)
A7BC3D13_cF192_219_de_b_c_ah-001.{ABCD}          cF192     F-43c               -         Mg3B7ClO13          #219      Boracite (part 2)   
AB11CD3_cP16_221_a_dg_b_c-001.{ABCD}             cP16      Pm-3m               -         CrFe11MoNi3         #221      Model of Ferrite (part 1)
AB27CD3_cP32_221_a_dij_b_c-001.{ABCD}            cP32      Pm-3m               -         CrFe2525Ni6         #221      Model of Austenite (part 1)
A5BCD6_cF416_228_eg_c_b_h-001.{ABCD}             cF416     Fd-3c               -         CuCrCl5[NH3]6       #228      CuCrCl56 (part 2)   
A2B3C12D3_cI160_230_a_c_h_d-001.{ABCD}           cI160     Ia-3d               S1_{4}    Co3Al2Si3O12        #230      Garnet (part 2)     

*****************************************************************************************************************************************************
ANRL - QUINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
A8B2C12D2E_oI50_23_bcfk_i_3k_j_a-001.{ABCDE}     oI50      I222                -         Cu8(Fe,Zn)3Sn2S12   #23       Stannoidite (part 2)
AB4C17D4E_tP54_90_a_g_c4g_g_c-001.{ABCDE}        tP54      P42_{1}2            -         BaCu4[VO][PO4]4     #90       BaCu44 (part 2)     
AB32CD4E8_tP184_93_i_16p_af_2p_4p-001.{ABCDE}    tP184     P4_{2}22            -         AsPh4CeS8P4Me8      #93       AsPh4CeS8P4Me8 (part 2)
A2B3C3DE7_hP48_145_2a_3a_3a_a_7a-001.{ABCDE}     hP48      P3_{2}              -         NaCa3[CO3]2F3[H2O]  #145      Sheldrickite (part 2)
A2BCD3E6_cF208_203_e_c_d_f_g-001.{ABCDE}         cF208     Fd-3                -         Na3Co(CO3)2Cl       #203      Pyrochlore (part 2) 
A4B2C6D16E_cF232_203_e_d_f_eg_a-001.{ABCDE}      cF232     Fd-3                -         Na6Mg2(SO4)(CO3)4   #203      Tychite (part 2)    
A6B2CD6E_cP64_208_m_ad_b_m_c-001.{ABCDE}         cP64      P4_{2}32            -         Cs2ZnFe[CN]6        #208      Cs2ZnFe6 (part 2)   

*****************************************************************************************************************************************************
ANRL - BCC-SQS

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB_aP16_2_4i_4i-001.AB                           aP16      P-1                 -         AB                  #2        bcc-SQS [10.1103/PhysRevB.69.214202] (O. Levy, part 3)
A5B11_mP16_6_2abc_2a3b3c-001.{AB}                mP16      Pm                  -         A5B11               #6        bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy, part 3)      
AB3_mC32_8_4a_12a-001.{AB}                       mC32      Cm                  -         AB3                 #8        bcc-SQS [10.1103/PhysRevB.69.214202] (O. Levy, part 3)
AB3_mC32_8_4a_4a4b-001.{AB}                      mC32      Cm                  -         AB3                 #8        bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy, part 3)
A3B13_oC32_38_ac_a2bcdef-001.{AB}                oC32      Amm2                -         A3B13               #38       bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy, part 3)
A3B5_oC32_38_abce_abcdf-001.{AB}                 oC32      Amm2                -         A3B5                #38       bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy, part 3)
AB7_hR16_166_c_c2h-001.{AB}                      hR16      R-3m                -         AB7                 #166      bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy, part 3)

*****************************************************************************************************************************************************
ANRL - miscellaneous

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2BCD4_tI16_82_ac_b_d_g-001.ABCD/BCAD/           tI16      I-4                 -         Cu2(Zn,Fe)SnS4      #82       Kesterite, 12 protos [Can. Mineral., 16, 131-137 (1978)] (C. Oses, part 3)
CABD/DBAC/BADC/ADBC/ACDB/CDAB/DACB/
DCBA/CBDA/BDCA   
A2B3_cI80_206_ad_e.{AB}                          tI80      Ia-3                -         NZn                 #206      NZn (ICSD #84918) (R. Friedrich)

*****************************************************************************************************************************************************
ANRL - fcc/sc derivative structures - BINARIES (Y. Lederer)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB3_mC8_12_a_di-001.{AB}                         mC8       C2/m                -         AB3                 #12       fcc derivative, equivalent to 16,18 [HTQC] 
AB_mC8_12_i_i-001.AB                             mC8       C2/m                -         AB                  #12       fcc derivative, equivalent to 17 [HTQC]
AB3_oP4_47_a_ct-001.{AB}                         oP4       Pmmm                -         AB3                 #47       fcc derivative, equivalent to 19,21 [HTQC]
A3B_oP4_47_cr_a-001.{AB}                         oP4       Pmmm                -         A3B                 #47       sc derivative
AB_oP4_51_e_e-001.AB                             oP4       Pmma                -         AB                  #51       sc derivative
AB_oP4_59_a_a-001.AB                             oP4       Pmmn                -         AB                  #59       fcc derivative, equivalent to 20 [HTQC]
AB2_oC6_65_a_i-001.{AB}                          oC6       Cmmm                -         AB2                 #65       sc derivative
A3B_oC8_65_ai_b-001.{AB}                         oC8       Cmmm                -         A3B                 #65       sc derivative
A3B_oC8_65_bj_a-001.{AB}                         oC8       Cmmm                -         A3B                 #65       sc derivative
AB_oC8_65_i_i-001.AB                             oC8       Cmmm                -         AB                  #65       sc derivative
AB_tP2_123_a_b-001.AB                            tP2       P4/mmm              -         AB                  #123      sc derivative
AB_tP2_123_a_c-001.AB                            tP2       P4/mmm              -         AB                  #123      sc derivative 
A2B_tP3_123_g_a-001.{AB}                         tP3       P4/mmm              -         A2B                 #123      sc derivative 
A3B_tP4_123_abc_d-001.{AB}                       tP4       P4/mmm              -         A3B                 #123      sc derivative 
A3B_tP4_123_ag_b-001.{AB}                        tP4       P4/mmm              -         A3B                 #123      sc derivative
A3B_tP4_123_cf_a-001.{AB}                        tP4       P4/mmm              -         A3B                 #123      sc derivative 
AB3_tP4_123_a_bh-001.{AB}                        tP4       P4/mmm              -         AB3                 #123      fcc derivative, equivalent to 13,15 [HTQC]
AB_tP4_123_g_g-001.AB                            tP4       P4/mmm              -         AB                  #123      sc derivative
A3B_tI8_139_ae_b-001.{AB}                        tI8       I4/mmm              -         A3B                 #139      sc derivative 
A2B_hP3_164_d_a-001.{AB}                         hP3       P-3m1               -         A2B                 #164      fcc derivative, equivalent to 9,10 [HTQC]
A2B_hP3_164_d_a-002.{AB}                         hP3       P-3m1               -         A2B                 #164      sc derivative 
A3B_hR4_166_bc_a-001.{AB}                        hR4       R-3m                -         A3B                 #166      sc derivative 
AB3_hR4_166_a_bc-001.{AB}                        hR4       R-3m                -         AB3                 #166      fcc derivative, equivalent to 27,29 [HTQC]
AB_hR4_166_c_c-001.AB                            hR4       R-3m                -         AB                  #166      fcc derivative, equivalent to 28 [HTQC]
AB_hR4_166_c_c-002.AB                            hR4       R-3m                -         AB                  #166      sc derivative

*****************************************************************************************************************************************************
ANRL - rocksalt/B2 derivative structures - TERNARIES (Y. Lederer)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB3C4_mC16_12_a_di_2i-001.{ABC}                  mC16      C2/m                -         AB3C4               #12       rocksalt derivative, AB3_mC8_12_a_di-001 with C "in between"
ABC2_mC16_12_i_i_adi-001.ABC/BCA/CAB             mC16      C2/m                -         ABC2                #12       rocksalt derivative, AB_mC8_12_i_i-001 with C "in between"
AB3C4_oP8_47_a_ct_egs-001.{ABC}                  oP8       Pmmm                -         AB3C4               #47       rocksalt derivative, AB3_oP4_47_a_ct-001 with C "in between"
A3BC4_oP8_47_eq_g_bdt-001.{ABC}                  oP8       Pmmm                -         A3BC4               #47       B2 derivative, A3B_oP4_47_cr_a-001 with C "in between", shifted origin
ABC2_oP8_51_e_e_2f-001.ABC/BCA/CAB               oP8       Pmma                -         ABC2                #51       B2 derivative, AB_oP4_51_e_e-001 with C "in between"
ABC2_oP8_59_a_a_2b-001.ABC/BCA/CAB               oP8       Pmmn                -         ABC2                #59       rocksalt derivative, AB_oP4_59_a_a-001 with C "in between"
ABC2_oC16_63_c_c_g-001.ABC/BCA/CAB               oC16      Cmcm                -         ABC2                #63       B2 derivative, AB_oC8_63_c_c-002 with C "in between" 
AB2C3_oC12_65_a_i_cj-001.{ABC}                   oC12      Cmmm                -         AB2C3               #65       B2 derivative, AB2_oC6_65_a_i-001 with C "in between" 
A3BC4_oC16_65_ai_b_q-001.{ABC}                   oC16      Cmmm                -         A3BC4               #65       B2 derivative, A3B_oC8_65_ai_b-001 with C "in between"
A3BC4_oC16_65_bj_a_eh-001.{ABC}                  oC16      Cmmm                -         A3BC4               #65       B2 derivative, A3B_oC8_65_bj_a-001 with C "in between"
AB3C4_oC16_65_a_bf_hi-001.{ABC}                  oC16      Cmmm                -         AB3C4               #65       rocksalt derivative, AB3_oC8_65_a_bf-001 with C "in between" 
ABC2_oC16_65_i_i_fh-001.ABC/BCA/CAB              oC16      Cmmm                -         ABC2                #65       B2 derivative, AB_oC8_65_i_i-001 with C "in between"
AB2C3_oI12_71_a_e_df-001.{ABC}                   oI12      Immm                -         AB2C3               #71       rocksalt derivative, AB2_oI6_71_a_g-001 with C "in between", shifted origin 
ABC2_tP4_123_a_b_h-001.ABC/BCA/CAB               tP4       P4/mmm              -         ABC2                #123      B2 derivative, AB_tP2_123_a_b with C "in between"; equivalent to TBCC011 [HTQC]
ABC2_tP4_123_a_c_e-001.ABC/BCA/CAB               tP4       P4/mmm              -         ABC2                #123      B2 derivative, AB_tP2_123_a_c-001 with C "in between"; equivalent to TBCC016 [HTQC]
ABC2_tP4_123_a_d_bc-001.ABC/BCA/CAB              tP4       P4/mmm              -         ABC2                #123      rocksalt derivative, AB_tP2_123_a_d-001 with C "in between" 
A2BC3_tP6_123_g_b_ch-001.{ABC}                   tP6       P4/mmm              -         A2BC3               #123      B2 derivative, A2B_tP3_123_g_a-001 with C "in between", shifted origin
A3BC4_tP8_123_abc_d_i-001.{ABC}                  tP8       P4/mmm              -         A3BC4               #123      B2 derivative, A3B_tP4_123_abc_d-001 with C "in between"
A3BC4_tP8_123_ag_b_2h-001.{ABC}                  tP8       P4/mmm              -         A3BC4               #123      B2 derivative, A3B_tP4_123_ag_b-001 with C "in between" 
A3BC4_tP8_123_cf_a_k-001.{ABC}                   tP8       P4/mmm              -         A3BC4               #123      B2 derivative, A3B_tP4_123_cf_a-001 with C "in between"
AB3C4_tP8_123_a_bh_cdg-001.{ABC}                 tP8       P4/mmm              -         AB3C4               #123      rocksalt derivative, AB3_tP4_123_a_bh-001 with C "in between"
ABC2_tP8_123_h_h_abg-001.ABC/BCA/CAB             tP8       P4/mmm              -         ABC2                #123      B2 derivative, AB_tP4_123_g_g-001 with C "in between", shifted origin
ABC2_tP8_129_c_c_2c-001.ABC/BCA/CAB              tP8       P4/nmm              -         ABC2                #129      rocksalt derivative, AB_tP4_129_c_c-002 with C "in between" 
AB2C3_tI12_139_a_e_be-001.{ABC}                  tI12      I4/mmm              -         AB2C3               #139      rocksalt derivative, AB_tI8_139_e_e-002 with C "in between"
A3BC4_tI16_139_ae_b_g-001.{ABC}                  tI16      I4/mmm              -         A3BC4               #139      B2 derivative, A3B_tI8_139_ae_b-001 with C "in between" 
AB3C4_tI16_139_a_bd_ce-001.{ABC}                 tI16      I4/mmm              -         AB3C4               #139      rocksalt derivative, A3B_tI8_139_bd_a-002.BA (swap A and B) with C "in between" 
ABC2_tI16_139_e_e_cd-001.ABC/BCA/CAB             tI16      I4/mmm              -         ABC2                #139      B2 derivative, AB_tI8_139_e_e-002 with C "in between" 
ABC2_tI16_141_a_b_e-001.ABC/BCA/CAB              tI16      I4_{1}/amd          -         ABC2                #141      rocksalt derivative, AB_tI8_141_a_b-002 with C "in between"
A2BC3_hP6_164_d_a_bd-001.{ABC}                   hP6       P-3m1               -         A2BC3               #164      rocksalt derivative, A2B_hP3_164_d_a-001 with C "in between"
A2BC3_hP6_164_d_b_ad-001.{ABC}                   hP6       P-3m1               -         A2BC3               #164      B2 derivative, A2B_hP3_164_d_a-002 with C "in between", shifted origin 
A3BC4_hR8_166_bc_a_2c-001.{ABC}                  hR8       R-3m                -         A3BC4               #166      B2 derivative, A3B_hR4_166_bc_a-001 with C "in between"
AB3C4_hR8_166_a_bc_2c-001.{ABC}                  hR8       R-3m                -         AB3C4               #166      rocksalt derivative, AB3_hR4_166_a_bc-001 with C "in between"
ABC2_hR8_166_c_c_abc-001.ABC/BCA/CAB             hR8       R-3m                -         ABC2                #166      rocksalt derivative, AB_hR4_166_c_c-001 with C "in between"
ABC2_hR8_166_c_c_abc-002.ABC/BCA/CAB             hR8       R-3m                -         ABC2                #166      B2 derivative, AB_hR4_166_c_c-002 with C "in between"
AB3C4_cP8_221_a_c_bd.{ABC}                       cP8       Pm-3m               -         AB3C4               #221      rocksalt derivative, AB3_cP4_221_a_c with C "in between" 

*****************************************************************************************************************************************************
ANRL - oxides - BINARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB_mC20_12_b2i_c2i-001.{AB}                      mC20      C2/m                -         OTi                 #12       OTi, ICSD #56694
A2B3_mC20_12_2i_3i-001.{AB}                      mC20      C2/m                -         Ga2O3               #12       Ga2O3, ICSD #83645, equivalent to beta-Ga2O3 (http://aflow.org/prototype-encyclopedia/A2B3_mC20_12_2i_3i.html, part 3)
A5B3_mC32_12_5i_3i-001.{AB}                      mC32      C2/m                -         O5Ti3               #12       O5Ti3, ICSD #75193
A3B_mP16_14_3e_e-001.{AB}                        mP16      P2_{1}/c            -         O3W                 #14       O3W, ICSD #84848
A2B3_mP20_14_2e_3e-001.{AB}                      mP20      P2_{1}/c            -         Bi2O3               #14       Bi2O3, ICSD #168806
AB_mC8_15_a_e-001.{AB}                           mC8       C2/c                -         CuO                 #15       CuO, ICSD #92368
A4B_mC20_15_2f_e-001.{AB}                        mC20      C2/c                -         O4Os                #15       O4Os, ICSD #23803
A2B5_oP28_19_2a_5a-001.{AB}                      oP28      P2_{1}2_{1}2_{1}    -         As2O5               #19       As2O5, ICSD #654040
A2B_oP24_33_4a_2a-001.{AB}                       oP24      Pna2_{1}            -         O2Sb                #33       O2Sb, ICSD #4109, equivalent to alpha-Sb2O4 (cervantite, http://aflow.org/prototype-encyclopedia/A2B_oP24_33_4a_2a.html, part 3)
AB3_oC16_40_b_a2b-001.{AB}                       oC16      Ama2                -         CrO3                #40       CrO3, ICSD #16031
A5B2_oP14_59_a2e_e-001.{AB}                      oP14      Pmmn                -         O5V2                #59       O5V2, ICSD #60767, equivalent to O5V2 (shcherbinaite (revised), http://aflow.org/prototype-encyclopedia/A5B2_oP14_59_a2f_f.html, part 3)
AB_oC16_64_e_f-001.{AB}                          oC16      Cmce                -         KO                  #64       KO, ICSD #180559
A4B3_tP28_135_gh_dh-001.{AB}                     tP28      P4_{2}/mbc          -         O4Pb3               #135      O4Pb3, ICSD #29094
A2B3_hP15_144_2a_3a-001.{AB}                     hP15      P3_{1}              -         B2O3                #144      B2O3, ICSD #24649
A2B_hR3_166_c_a-001.{AB}                         hR3       R-3m                -         Cs2O                #166      Cs2O, ICSD #27919
AB2_hR6_166_c_2c-001.{AB}                        hR6       R-3m                -         OTl2                #166      OTl2, ICSD #16220
AB_hP12_189_fg_eh-001.{AB}                       hP12      P-62m               -         NaO                 #189      NaO, ICSD #25526
AB_hP8_194_ac_f-001.{AB}                         hP8       P6_{3}/mmc          -         LiO                 #194      LiO, ICSD #152183
A2B3_cI80_199_a2b_2c-001.{AB}                    cI80      I2_{1}3             -         In2O3               #199      In2O3, ICSD #33649
A3B2_cF80_227_f_e-001.{AB}                       cF80      Fd-3m               -         O3Sb2               #227      O3Sb2, ICSD #31102

*****************************************************************************************************************************************************
ANRL - oxides - TERNARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A4B4C_aP18_2_4i_4i_i-001.{ABC}                   aP18      P-1                 -         Na4O4Si             #2        Na4O4Si, ICSD #15500
A2B7C2_aP22_2_2i_7i_2i-001.{ABC}                 aP22      P-1                 -         Ca2O7V2             #2        Ca2O7V2, ICSD #421266
A2B7C2_aP22_2_2i_7i_2i-002.{ABC}                 aP22      P-1                 -         Mg2O7V2             #2        Mg2O7V2, ICSD #2321
AB3C_aP30_2_3i_9i_3i-001.{ABC}                   aP30      P-1                 -         CaO3Si              #2        CaO3Si, ICSD #201537, equivalent to CaO3Si (wollastonite, http://aflow.org/prototype-encyclopedia/AB3C_aP30_2_3i_9i_3i.html, part 3)
A2B5C_aP32_2_4i_10i_2i-001.{ABC}                 aP32      P-1                 S0_{1}    Al2O5Si             #2        Al2O5Si, ICSD #85742, equivalent to Al2O5Si (kyanite, http://aflow.org/prototype-encyclopedia/A2B5C_aP32_2_4i_10i_2i.html, part 3)
A2B4C_mP28_4_4a_8a_2a-001.{ABC}                  mP28      P2_{1}              -         Al2O4Sr             #4        Al2O4Sr, ICSD #160297
A3BC_mC60_5_ab8c_ab2c_3c-001.{ABC}               mC60      C2                  -         O3SiSr              #5        O3SiSr, ICSD #32542
A3B5C_mC54_8_3a3b_9a3b_3a-001.{ABC}              mC54      Cm                  -         Ca3O5Si             #8        Ca3O5Si, ICSD #81100
A2B7C3_mP24_11_2e_7e_3e-001.{ABC}                mP24      P2_{1}/m            -         Na2O7Ti3            #11       Na2O7Ti3, ICSD #15463
AB6C2_mC18_12_a_3i_i-001.{ABC}                   mC18      C2/m                -         CaO6V2              #12       CaO6V2, ICSD #166516
AB6C2_mC18_12_a_3i_i-002.{ABC}                   mC18      C2/m                -         MgO6V2              #12       MgO6V2, ICSD #10391
ABC4_mC48_12_gi_hi_2i3j-001.{ABC}                mC48      C2/m                -         MgMoO4              #12       MgMoO4, ICSD #20418
ABC4_mC48_12_gi_hi_2i3j-002.{ABC}                mC48      C2/m                -         FeMoO4;             #12       FeMoO4, ICSD #43012
AB4C_mP12_13_f_2g_e-001.{ABC}                    mP12      P2/c                H0_{6}    MgO4W               #13       MgO4W, ICSD #67903, equivalent to MgO4W (huanzalaite, H0_{6}, http://aflow.org/prototype-encyclopedia/AB4C_mP12_13_f_2g_e.html, part 3)
AB4C_mP12_13_f_2g_e-002.{ABC}                    mP12      P2/c                -         NiO4W               #13       NiO4W, ICSD #16685
AB3C_mP40_14_2e_6e_2e-001.{ABC}                  mP40      P2_{1}/c            -         MgO3Si              #14       MgO3Si, ICSD #30895
A2B3C_mC24_15_2e_af_e-001.{ABC}                  mC24      C2/c                -         Li2O3Zr             #15       Li2O3Zr, ICSD #94893
AB3C_mC40_15_2e_3f_f-001.{ABC}                   mC40      C2/c                -         NaO3V               #15       NaO3V, ICSD #2103
A2B3C_mC48_15_aef_3f_2e-001.{ABC}                mC48      C2/c                -         Li2O3Ti             #15       Li2O3Ti, ICSD #162215, equivalent to Na2O3Pr (http://aflow.org/prototype-encyclopedia/A2B3C_mC48_15_aef_3f_2e.html, part 3)
A4BC7_mC48_15_2f_e_e3f-001.{ABC}                 mC48      C2/c                -         Al4CaO7             #15       Al4CaO7, ICSD #14270
AB3C_mC60_15_cf_e4f_ef-001.{ABC}                 mC60      C2/c                -         CaO3Si              #15       CaO3Si, ICSD #87694
ABC2_oP16_33_a_a_2a-001.{ABC}                    oP16      Pna2_{1}            -         FeNaO2;             #33       FeNaO2, ICSD #186309
A2B3C_oC24_36_b_ab_a-001.{ABC}                   oC24      Cmc2_{1}            -         Li2O3Si             #36       Li2O3Si, ICSD #28192
A2B3C_oC24_36_b_ab_a-002.{ABC}                   oC24      Cmc2_{1}            -         Na2O3Si             #36       Na2O3Si, ICSD #24664
A2B5C_oP32_58_eg_3gh_g-001.{ABC}                 oP32      Pnnm                S0_{2}    Al2O5Si             #58       Al2O5Si, ICSD #30679, equivalent to Al2O5Si (http://aflow.org/prototype-encyclopedia/A2B5C_oP32_58_eg_3gh_g.html, part 3)
AB2C4_oC28_63_c_ac_fg-001.{ABC}                  oC28      Cmcm                H1_{8}    CrNa2O4             #63       CrNa2O4, ICSD #76001, equivalent to CrNa2O4 (http://aflow.org/prototype-encyclopedia/AB2C4_oC28_63_c_bc_fg.html, part 3)
AB5C2_oC32_63_c_c2f_f-001.{ABC}                  oC32      Cmcm                E4_{1}    MgO5Ti2/Fe2O5Ti     #63       MgO5Ti2, ICSD #37232, equivalent to Fe2O5Ti (pseudobrookite, http://aflow.org/prototype-encyclopedia/A2B5C_oC32_63_f_c2f_c.html, part 3)
ABC4_tI24_88_a_b_f-001.{ABC}                     tI24      I4_{1}/a            -         BaMoO4              #88       BaMoO4, ICSD #50821
ABC4_tI24_88_a_b_f-002.{ABC}                     tI24      I4_{1}/a            H0_{4}    CaMoO4/CaO4W        #88       CaMoO4, ICSD #417513, equivalent to CaO4W (scheelite, H0_{4}, http://aflow.org/prototype-encyclopedia/AB4C_tI24_88_b_f_a.html, part 3)
A4BC2_tP28_91_2d_b_ac-001.{ABC}                  tP28      P4_{1}22            -         O4TiZn2             #91       O4TiZn2, ICSD #109093
A2BC4_hP56_173_2b2c_ac_b5c-001.{ABC}             hP56      P6_{3}              -         Al2BaO4             #173      Al2BaO4, ICSD #246027
AB2C4_cF56_227_b_c_e-001.{ABC}                   cF56      Fd-3m               -         MoNa2O4             #227      MoNa2O4, ICSD #151970
    
*****************************************************************************************************************************************************
ANRL - nitrides - BINARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2B3_cI80_206_ad_e-001.{AB}                      cI80      Ia-3                -         N2Zn3               #206      N2Zn3, ICSD #84918

*****************************************************************************************************************************************************
ANRL - nitrides - TERNARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A3B3C_mC56_9_6a_6a_2a-001.{ABC}                  mC56      Cc                  -         N3Na3W1             #9        N3Na3W1, ICSD #75364
ABC_oP12_62_c_c_c-001.{ABC}                      oP12      Pnma                -         CaLiN               #62       CaLiN, ICSD #107304
A2B2C_tI10_139_e_e_a-001.{ABC}                   tI10      I4/mmm              -         Ca2N2Zn             #139      Ca2N2Zn, ICSD #69049
A3BC2_hP6_191_f_a_d-001.{ABC}                    hP6       P6/mmm              -         Na3TaTi2            #191      Na3TaTi2, ICSD #186418

*****************************************************************************************************************************************************
ANRL - metal nitrides - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2B_mP6_6_3ab_ab-001.{AB}                        mP6       Pm                  -         N2Re1               #6        metal-nitride; N2Re1, ICSD #187448
A2B_mC12_8_4a_2a-001.{AB}                        mC12      Cm                  -         N2Re1               #8        metal-nitride; N2Re1, ICSD #187447
A3B_mC16_8_6a_b-001.{AB}                         mC16      Cm                  -         N3Rb1               #8        metal-nitride; N3Rb1, ICSD #155169
A2B_mP6_11_2e_e-001.{AB}                         mP6       P2_{1}/m            -         N2Re1               #11       metal-nitride; N2Re1, ICSD #187444
AB6_mP14_11_e_6e-001.{AB}                        mP14      P2_{1}/m            -         Ba1N6               #11       metal-nitride; Ba1N6, ICSD #14244
A3B_mC8_12_di_a-001.{AB}                         mC8       C2/m                -         N3Na1               #12       metal-nitride; N3Na1, ICSD #29370
A3B_mC8_12_di_a-002.{AB}                         mC8       C2/m                -         N3Na1               #12       metal-nitride; N3Na1, ICSD #29376
AB3_mC8_12_a_di-002.{AB}                         mC8       C2/m                -         Li1N3               #12       metal-nitride; Li1N3, ICSD #34675
AB3_mC8_12_a_di-003.{AB}                         mC8       C2/m                -         Li1N3               #12       metal-nitride; Li1N3, ICSD #181559
A2B_mC12_12_2i_i-001.{AB}                        mC12      C2/m                -         N2Re1               #12       metal-nitride; N2Re1, ICSD #187441
AB_mC16_12_2i_2i-001.{AB}                        mC16      C2/m                -         N1Sr1               #12       metal-nitride; N1Sr1, ICSD #411555
A3B2_mC30_12_a4i_3i-001.{AB}                     mC30      C2/m                -         Ca3N2               #12       metal-nitride; Ca3N2, ICSD #162794
A3B2_mC30_12_a4i_3i-002.{AB}                     mC30      C2/m                -         Ca3N2               #12       metal-nitride; Ca3N2, ICSD #169726
A2B_mP6_13_g_e-001.{AB}                          mP6       P2/c                -         N2Re1               #13       metal-nitride; N2Re1, ICSD #187451
AB2_mP12_14_e_2e-001.{AB}                        mP12      P2_{1}/c            -         Ir1N2               #14       metal-nitride; Ir1N2, ICSD #160623
A2B_mP12_14_2e_e-004.{AB}                        mP12      P2_{1}/c            -         N2Os1               #14       metal-nitride; N2Os1, ICSD #240759
AB3_mP16_14_e_3e-001.{AB}                        mP16      P2_{1}/c            -         Hg2N6               #14       metal-nitride; Hg2N6, ICSD #98661
AB2_mC12_15_e_f-001.{AB}                         mC12      C2/c                C_{g}     Ba1N2               #15       metal-nitride; Ba1N2, ICSD #280681, equivalent to C2Th (C_{g}, http://aflow.org/prototype-encyclopedia/A2B_mC12_15_f_e.html, part 3)
AB_oC16_36_2a_2a-001.{AB}                        oC16      Cmc2_{1}            -         N1Os1               #36       metal-nitride; N1Os1, ICSD #167514
AB_oP8_53_h_h-001.{AB}                           oP8       Pmna                -         N1Os1               #53       metal-nitride; N1Os1, ICSD #167513
A2B_oP12_57_e_c-001.{AB}                         oP12      Pbcm                -         N2Re1               #57       metal-nitride; N2Re1, ICSD #187450
A2B_oP6_58_g_a-002.{AB}                          oP6       Pnnm                -         N2Os1               #58       metal-nitride; N2Os1, ICSD #157283
A2B_oP6_58_g_a-003.{AB}                          oP6       Pnnm                -         N2Pt1               #58       metal-nitride; N2Pt1, ICSD #166463
AB2_oP6_58_a_g-004.{AB}                          oP6       Pnnm                -         Ir1N2               #58       metal-nitride; Ir1N2, ICSD #160620
AB2_oP6_58_a_g-005.{AB}                          oP6       Pnnm                -         Au1N2               #58       metal-nitride; Au1N2, ICSD #166465
A2B_oP12_60_d_c-001.{AB}                         oP12      Pbcn                -         Fe2N1               #60       metal-nitride; Fe2N1, ICSD #150889 similar to 1213.AB,1213.BA,1222.AB,1222.BA,1284.AB,1284.BA [HTQC], equivalent to zeta-Fe2N (http://aflow.org/prototype-encyclopedia/A2B_oP12_60_d_c.Fe2N.html, part 3)
A2B_oP12_62_2c_c-005.{AB}                        oP12      Pnma                -         N2Re1               #62       metal-nitride; N2Re1, ICSD #187446
AB3_oP16_62_c_cd-002.{AB}                        oP16      Pnma                -         N1Na3               #62       metal-nitride; N1Na3, ICSD #165990
A3B_oP16_62_cd_c-002.{AB}                        oP16      Pnma                -         Fe3N1               #62       metal-nitride; Fe3N1, ICSD #260758 similar to 1076.AB,1076.BA [HTQC]
AB_oC8_63_c_c-003.AB                             oC8       Cmcm                -         Al1N1               #63       metal-nitride; Al1N1, ICSD #163951
A5B3_oC32_63_c2f_cf-001.{AB}                     oC32      Cmcm                -         N5Ta3               #63       metal-nitride; N5Ta3, ICSD #16253
A5B3_oC32_63_c2f_cf-002.{AB}                     oC32      Cmcm                -         N5Ta3               #63       metal-nitride; N5Ta3, ICSD #66533
A2B_oC12_65_gj_e-001.{AB}                        oC12      Cmmm                -         N2Re1               #65       metal-nitride; N2Re1, ICSD #187452
AB_oI4_71_b_c-001.AB                             oI4       Immm                -         Cr1N1               #71       metal-nitride; Cr1N1, ICSD #53146
A2B_oI12_71_2e_f-001.{AB}                        oI12      Immm                -         N2Re1               #71       metal-nitride; N2Re1, ICSD #187449
AB_oI16_71_abe_n-001.{AB}                        oI16      Immm                -         Li2N2               #71       metal-nitride; Li2N2, ICSD #423831
AB3_oI16_72_b_cj-001.{AB}                        oI16      Ibam                -         Ag1N3               #72       metal-nitride; Ag1N3, ICSD #27135
A5B4_tI18_87_ah_h-001.{AB}                       tI18      I4/m                -         N5Nb4               #87       metal-nitride; N5Nb4, ICSD #26251
AB3_tI32_88_d_cf-001.{AB}                        tI32      I4_{1}/a            -         Cu1N3               #88       metal-nitride; Cu1N3, ICSD #30633, equivalent to Copper (I) Azide (CuN3, http://aflow.org/prototype-encyclopedia/AB3_tI32_88_d_cf.html, part 3)
AB_tP2_123_d_a-001.AB                            tP2       P4/mmm              -         N1Pr1               #123      metal-nitride; N1Pr1, ICSD #168645
A3B_tP4_123_ag_d-001.{AB}                        tP4       P4/mmm              -         N3Rb1               #123      metal-nitride; N3Rb1, ICSD #16963
A3B_tP8_123_egh_ab-001.{AB}                      tP8       P4/mmm              -         Cu3N1               #123      metal-nitride; Cu3N1, ICSD #180237
A2B_tP6_127_g_b-001.{AB}                         tP6       P4/mbm              -         N2Re1               #127      metal-nitride; N2Re1, ICSD #187442
AB_tP4_131_e_c-001.{AB}                          tP4       P4_{2}/mmc          -         N1Tc1               #131      metal-nitride; N1Tc1, ICSD #187709
AB_tI4_139_a_b-001.AB                            tI4       I4/mmm              -         Mn1N1               #139      metal-nitride; Mn1N1, ICSD #106932
A3B2_tI10_139_ae_e-001.{AB}                      tI10      I4/mmm              -         Mn3N2               #139      metal-nitride; Mn3N2, ICSD #71638
A3B2_tI10_139_ae_e-002.{AB}                      tI10      I4/mmm              -         Mn3N2               #139      metal-nitride; Mn3N2, ICSD #84202
A3B4_tI14_139_ad_ce-001.{AB}                     tI14      I4/mmm              -         N3Nb4               #139      metal-nitride; N3Nb4, ICSD #76389
A8B_tI18_139_deh_a-001.{AB}                      tI18      I4/mmm              D2_{g}    Fe16N2              #139      metal-nitride; Fe16N2, ICSD #41953, equivalent to Fe8N (D2_{g}, http://aflow.org/prototype-encyclopedia/A8B_tI18_139_deh_a.html, part 3)
AB3_tI16_140_a_dh-001.{AB}                       tI16      I4/mcm              -         K1N3                #140      metal-nitride; K1N3, ICSD #24007
AB3_tI16_140_a_dh-002.{AB}                       tI16      I4/mcm              -         Cs1N3               #140      metal-nitride; Cs1N3, ICSD #25008
AB3_tI16_140_a_dh-003.{AB}                       tI16      I4/mcm              -         Ag1N3               #140      metal-nitride; Ag1N3, ICSD #183201
AB3_tI16_140_a_ch-001.{AB}                       tI16      I4/mcm              -         Cu1N3               #140      metal-nitride; Cu1N3, ICSD #187018
AB3_tI16_140_c_ah-001.{AB}                       tI16      I4/mcm              -         K1N3                #140      metal-nitride; K1N3, ICSD #1145
AB2_tI12_141_a_e-001.{AB}                        tI12      I4_{1}/amd          -         N1Ti2               #141      metal-nitride; N1Ti2, ICSD #23403
A2B_tI12_141_e_a-003.{AB}                        tI12      I4_{1}/amd          -         Mo2N1               #141      metal-nitride; Mo2N1, ICSD #30593 similar to 725.AB,725.BA [HTQC]
A3B_hR4_160_3a_a-001.{AB}                        hR4       R3m                 -         N3Na1               #160      metal-nitride; N3Na1, ICSD #644523
A3B2_hR5_160_3a_2a-001.{AB}                      hR5       R3m                 -         Be3N2               #160      metal-nitride; Be3N2, ICSD #185490
AB_hP16_162_ek_ci-001.{AB}                       hP16      P-31m               -         Mo1N1               #162      metal-nitride; Mo1N1, ICSD #43559
A2B_hP3_164_d_a-003.{AB}                         hP3       P-3m1               -         N2Zr1               #164      metal-nitride; N2Zr1, ICSD #262746 similar to 190,191 [HTQC]
A3B2_hP5_164_ad_d-002.{AB}                       hP5       P-3m1               -         Ca3N2               #164      metal-nitride; Ca3N2, ICSD #162795
A3B2_hP5_164_ad_d-003.{AB}                       hP5       P-3m1               -         Ca3N2               #164      metal-nitride; Ca3N2, ICSD #162796
A3B2_hP5_164_ad_d-004.{AB}                       hP5       P-3m1               -         Ca3N2               #164      metal-nitride; Ca3N2, ICSD #169727
A3B2_hP5_164_ac_d-001.{AB}                       hP5       P-3m1               -         N3V2                #164      metal-nitride; N3V2, ICSD #182700
AB_hP16_164_ci_di-001.{AB}                       hP16      P-3m1               -         Mo1N1               #164      metal-nitride; Mo1N1, ICSD #60168
A2B_hR3_166_c_a-002.{AB}                         hR3       R-3m                -         Ca2N1               #166      metal-nitride; Ca2N1, ICSD #22231
A2B_hR3_166_c_a-003.{AB}                         hR3       R-3m                -         Ba2N1               #166      metal-nitride; Ba2N1, ICSD #409851
AB2_hR3_166_a_c-001.{AB}                         hR3       R-3m                -         N1Sr2               #166      metal-nitride; N1Sr2, ICSD #23530 similar to 692.AB,692.BA [HTQC]
A3B_hR4_166_bc_a-002.{AB}                        hR4       R-3m                -         N3Na1               #166      metal-nitride; N3Na1, ICSD #1144
AB_hP4_186_b_b-002.AB                            hP4       P6_{3}mc            -         Ga1N1               #186      metal-nitride; Ga1N1, ICSD #153888
AB_hP4_186_b_b-003.AB                            hP4       P6_{3}mc            -         La1N1               #186      metal-nitride; La1N1, ICSD #162195
AB_hP4_186_b_b-004.AB                            hP4       P6_{3}mc            -         Mo1N1               #186      metal-nitride; Mo1N1, ICSD #168369
AB_hP4_186_b_a-002.{AB}                          hP4       P6_{3}mc            -         Ga1N1               #186      metal-nitride; Ga1N1, ICSD #159250
AB_hP16_186_ac_bc-001.{AB}                       hP16      P6_{3}mc            -         Mo1N1               #186      metal-nitride; Mo1N1, ICSD #76280
AB_hP2_187_d_a-002.AB                            hP2       P-6m2               -         Ir1N1               #187      metal-nitride; Ir1N1, ICSD #186245 similar to 1015.AB [HTQC]
AB2_hP3_187_a_h-001.{AB}                         hP3       P-6m2               -         Hf1N2               #187      metal-nitride; Hf1N2, ICSD #290427
A2B_hP3_187_h_a-001.{AB}                         hP3       P-6m2               -         N2W1                #187      metal-nitride; N2W1, ICSD #290433 similar to 695.AB,695.BA,1219.AB,1219.BA [HTQC]
AB_hP6_189_f_ad-001.{AB}                         hP6       P-62m               -         N1Ta1               #189      metal-nitride; N1Ta1, ICSD #1396 similar to 721.AB,721.BA [HTQC]
A2B_hP3_191_e_a-001.{AB}                         hP3       P6/mmm              -         N2Os1               #191      metal-nitride; N2Os1, ICSD #260545
AB_hP6_191_f_ad-002.{AB}                         hP6       P6/mmm              -         N1Ta1               #191      metal-nitride; N1Ta1, ICSD #25659
A2B_hP6_191_af_d-001.{AB}                        hP6       P6/mmm              -         N4Ta2               #191      metal-nitride; N4Ta2, ICSD #182351
AB_hP6_191_bc_ad-001.AB                          hP6       P6/mmm              -         N1Y1                #191      metal-nitride; N1Y1, ICSD #161079
A3B_hP8_193_g_b-001.{AB}                         hP8       P6_{3}/mcm          -         Ba3N1               #193      metal-nitride; Ba3N1, ICSD #77730
AB_hP4_194_c_a-002.{AB}                          hP4       P6_{3}/mmc          -         N1Nb1               #194      metal-nitride; N1Nb1, ICSD #76384 similar to 1031.AB,1031.BA,1232.AB,1232.BA [HTQC]
AB_hP4_194_c_a-003.{AB}                          hP4       P6_{3}/mmc          -/L'3_{0} N1Ta1               #194      metal-nitride; N1Ta1, ICSD #105123 similar to 1231.AB,1231.BA [HTQC], equivalent to Fe2N (L'3_{0}, http://aflow.org/prototype-encyclopedia/AB_hP4_194_c_a.Fe2N.html, part 3)
AB_hP4_194_c_a-004.{AB}                          hP4       P6_{3}/mmc          -         N1Re1               #194      metal-nitride; N1Re1, ICSD #181300 similar to 1030.AB,1030.BA [HTQC]
AB_hP4_194_c_a-005.{AB}                          hP4       P6_{3}/mmc          -         N1Nb1               #194      metal-nitride; N1Nb1, ICSD #185559 similar to 1030.AB,1030.BA,1232.AB,1232.BA [HTQC]
AB_hP4_194_c_d-002.AB                            hP4       P6_{3}/mmc          -         Al1N1               #194      metal-nitride; Al1N1, ICSD #163950
AB_hP4_194_a_c-001.{AB}                          hP4       P6_{3}/mmc          -         Cd1N1               #194      metal-nitride; Cd1N1, ICSD #185566 similar to 1031.AB,1031.BA,1231.AB,1231.BA [HTQC]
AB2_hP6_194_c_f-002.{AB}                         hP6       P6_{3}/mmc          -         N1Re2               #194      metal-nitride; N1Re2, ICSD #169882 similar to 1243.AB,1243.BA [HTQC]
AB2_hP6_194_c_e-001.{AB}                         hP6       P6_{3}/mmc          -         Hf1N2               #194      metal-nitride; Hf1N2, ICSD #290428
A2B_hP6_194_e_c-001.{AB}                         hP6       P6_{3}/mmc          -         N2Ta1               #194      metal-nitride; N2Ta1, ICSD #290431 similar to 696.AB,696.BA [HTQC]
AB_hP8_194_ac_f-002.{AB}                         hP8       P6_{3}/mmc          -         N1Nb1               #194      metal-nitride; N1Nb1, ICSD #76008 similar to 1264.AB,1264.BA [HTQC]
A3B2_hP10_194_bf_ac-001.{AB}                     hP10      P6_{3}/mmc          -         Be3N2               #194      metal-nitride; Be3N2, ICSD #25656
A3B2_hP10_194_bf_ac-002.{AB}                     hP10      P6_{3}/mmc          -         Ca3N2               #194      metal-nitride; Ca3N2, ICSD #162797
A3B2_hP10_194_cf_f-001.{AB}                      hP10      P6_{3}/mmc          -         N3W2                #194      metal-nitride; N3W2, ICSD #186207
AB_hP16_194_bh_ag-001.{AB}                       hP16      P6_{3}/mmc          -         Mo1N1               #194      metal-nitride; Mo1N1, ICSD #106926
AB3_cI32_204_c_g-001.{AB}                        cI32      Im-3                -         Co1N3               #204      metal-nitride; Co1N3, ICSD #162105
AB3_cI32_204_c_g-002.{AB}                        cI32      Im-3                -         Co1N3               #204      metal-nitride; Co1N3, ICSD #162106
A4B_cP5_221_ac_b.{AB}                            cP5       Pm-3m               -         Mn4N1               #221      metal-nitride; Mn4N1, ICSD #44369 similar to 679.AB,679.BA [HTQC]
AB4_cP5_221_a_bc.{AB}                            cP5       Pm-3m               L'1_{0}   N1Ni4 (Fe4N)        #221      metal-nitride; N1Ni4, ICSD #76403, equivalent to Fe4N (http://aflow.org/prototype-encyclopedia/A4B_cP5_221_bc_a.html, part 3)
A4B3_cP7_221_ac_d.{AB}                           cP7       Pm-3m               -         N4W3                #221      metal-nitride; N4W3, ICSD #30370 similar to 708.AB,708.BA,1256.AB,1256.BA [HTQC]

*****************************************************************************************************************************************************
ANRL - metal nitrides - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB2C5_aP8_1_a_2a_5a-001.{ABC}                    aP8       P1                  -         Li1N2Na5            #1        metal-nitride; Li1N2Na5, ICSD #92316
AB2C5_mC16_5_a_c_a2bc-001.{ABC}                  mC16      C2                  -         Li1N2Na5            #5        metal-nitride; Li1N2Na5, ICSD #92315
A3B2C3_mP8_6_2ab_ab_a2b-001.{ABC}                mP8       Pm                  -         Li3N2Na3            #6        metal-nitride; Li3N2Na3, ICSD #92312
A5B2C_mP8_6_3a2b_ab_b-001.{ABC}                  mP8       Pm                  -         Li5N2Na1            #6        metal-nitride; Li5N2Na1, ICSD #92314
A2BC2_mC20_12_2i_i_2i-001.{ABC}                  mC20      C2/m                -         Ca2Fe1N2            #12       metal-nitride; Ca2Fe1N2, ICSD #72389
AB2C12_mC30_12_a_i_6i-001.{ABC}                  mC30      C2/m                -         Cd1K2N12            #12       metal-nitride; Cd1K2N12, ICSD #31297
AB2C12_mC30_12_a_i_6i-002.{ABC}                  mC30      C2/m                -         Cd1K2N12            #12       metal-nitride; Cd1K2N12, ICSD #659621
ABC_mP12_14_e_e_e-001.{ABC}                      mP12      P2_{1}/c            -         Be1Li1N1            #14       metal-nitride; Be1Li1N1, ICSD #402341
ABC_oP12_62_c_c_c-002.{ABC}                      oP12      Pnma                -         Li1Mg1N1            #62       metal-nitride; Li1Mg1N1, ICSD #93259
AB2C_oC16_63_c_2c_c-002.{ABC}                    oC16      Cmcm                -         Au1Ca2N1            #63       metal-nitride; Au1Ca2N1, ICSD #85528
A2BC_oC16_63_2c_c_c-001.{ABC}                    oC16      Cmcm                -         Ca2In1N1            #63       metal-nitride; Ca2In1N1, ICSD #96228
ABC3_oC20_63_c_a_cf-001.{ABC}                    oC20      Cmcm                -         Al1N1Zr3            #63       metal-nitride; Al1N1Zr3, ICSD #29521
A3BC3_oC28_63_cg_c_cg-001.{ABC}                  oC28      Cmcm                -         Ca3Cr1N3            #63       metal-nitride; Ca3Cr1N3, ICSD #40205
A2BC_oC8_65_i_c_b-001.{ABC}                      oC8       Cmmm                -         Li4N2Na2            #65       metal-nitride; Li4N2Na2, ICSD #92307
A7B2C3_oC24_65_a2ij_i_cj-001.{ABC}               oC24      Cmmm                -         Ca7N2Tl3            #65       metal-nitride; Ca7N2Tl3, ICSD #98178
A2BC3_oI12_71_e_a_de-001.{ABC}                   oI12      Immm                -         Ce2Mn1N3            #71       metal-nitride; Ce2Mn1N3, ICSD #50579
ABC_tP6_99_ab_ab_ab-001.{ABC}                    tP6       P4mm                -         Cr1N1Nb1            #99       metal-nitride; Cr1N1Nb1, ICSD #23779
A3B2C3_tP8_115_abc_g_de-001.{ABC}                tP8       P-4m2               -         Li3N2Na3            #115      metal-nitride; Li3N2Na3, ICSD #92311
ABC_tP3_123_c_a_b-001.{ABC}                      tP3       P4/mmm              -         Fe1N1Ni1            #123      metal-nitride; Fe1N1Ni1, ICSD #53505
A5B2C_tP8_123_ai_bc_d-001.{ABC}                  tP8       P4/mmm              -         Li5N2Na1            #123      metal-nitride; Li5N2Na1, ICSD #92313
ABC_tP6_129_c_c_c-001.{ABC}                      tP6       P4/nmm              -         Ca1Ga1N1            #129      metal-nitride; Ca1Ga1N1, ICSD #2027
ABC_tP6_129_c_c_c-002.{ABC}                      tP6       P4/nmm              -         Mo1N1Ta1            #129      metal-nitride; Mo1N1Ta1, ICSD #100437
ABC2_tP8_129_c_c_bc-001.{ABC}                    tP8       P4/nmm              -         Ba1Hf1N2            #129      metal-nitride; Ba1Hf1N2, ICSD #50994
ABC2_tP8_129_a_c_bc-001.{ABC}                    tP8       P4/nmm              -         Li2N2Na4            #129      metal-nitride; Li2N2Na4, ICSD #92309
ABC_tP6_131_e_c_b-001.{ABC}                      tP6       P4_{2}/mmc          -         Ca1N1Ni1            #131      metal-nitride; Ca1N1Ni1, ICSD #69044
ABC_tP6_131_b_c_e-001.{ABC}                      tP6       P4_{2}/mmc          -         Li1N1Sr1            #131      metal-nitride; Li1N1Sr1, ICSD #87414
A2BC_tI8_139_d_b_a-001.ABC/BCA/CAB               tI8       I4/mmm              -         Li4N2Na2            #139      metal-nitride; Li4N2Na2, ICSD #92305
A2B2C_tI10_139_e_e_a-002.{ABC}                   tI10      I4/mmm              -         N2Sr2Zn1            #139      metal-nitride; N2Sr2Zn1, ICSD #80376
A2B2C_tI10_139_e_e_a-003.{ABC}                   tI10      I4/mmm              -         Ba2N2Zn1            #139      metal-nitride; Ba2N2Zn1, ICSD #80377
A3B6C2_tI22_139_ae_eg_e-001.{ABC}                tI22      I4/mmm              -         La3N6V2             #139      metal-nitride; La3N6V2, ICSD #98477
A3B6C2_tI22_139_ae_eg_e-002.{ABC}                tI22      I4/mmm              -         La3N6Nb2            #139      metal-nitride; La3N6Nb2, ICSD #411473
A2B2C_tI20_140_h_h_a-001.{ABC}                   tI20      I4/mcm              -         Be4N4Sr2            #140      metal-nitride; Be4N4Sr2, ICSD #413356
A2BC2_tI20_140_h_a_h-001.{ABC}                   tI20      I4/mcm              -         Be4Ca2N4            #140      metal-nitride; Be4Ca2N4, ICSD #413357
AB2C2_tI20_140_a_h_h-001.{ABC}                   tI20      I4/mcm              -         Ba1Be2N2            #140      metal-nitride; Ba1Be2N2, ICSD #415304
AB2C_hR4_160_a_2a_a-001.{ABC}                    hR4       R3m                 -         Cr1N2W1             #160      metal-nitride; Cr1N2W1, ICSD #84639
A2B2C_hP5_164_d_d_a-001.{ABC}                    hP5       P-3m1               -         Li2N2Zr1            #164      metal-nitride; Li2N2Zr1, ICSD #16231
A2B2C_hP5_164_d_d_a-002.{ABC}                    hP5       P-3m1               -         Mg2N2Sr1            #164      metal-nitride; Mg2N2Sr1, ICSD #410826, equivalent to Ce2O2S (http://aflow.org/prototype-encyclopedia/A2B2C_hP5_164_d_d_a.html, part 3)
AB2C2_hP5_164_a_d_d-001.{ABC}                    hP5       P-3m1               -         Ce1Li2N2            #164      metal-nitride; Ce1Li2N2, ICSD #34003
AB2C_hR4_166_a_c_b-001.{ABC}                     hR4       R-3m                -         Cu1N2Ta1            #166      metal-nitride; Cu1N2Ta1, ICSD #71136
A3BC3_hP14_176_h_c_h-001.{ABC}                   hP14      P6_{3}/m            -         Ba3Fe1N3            #176      metal-nitride; Ba3Fe1N3, ICSD #36502
ABC_hP6_186_b_b_a-001.{ABC}                      hP6       P6_{3}mc            -         N1Na1Sn1            #186      metal-nitride; N1Na1Sn1, ICSD #172471
AB2C_hP8_186_a_ab_b-001.{ABC}                    hP8       P6_{3}mc            -         Mn1N2W1             #186      metal-nitride; Mn1N2W1, ICSD #80029
ABC_hP3_187_a_c_d-001.{ABC}                      hP3       P-6m2               -         Li1N1Ni1            #187      metal-nitride; Li1N1Ni1, ICSD #247028
A5B3C3_hP11_189_dg_g_f-001.{ABC}                 hP11      P-62m               -         Li5N3Ni3            #189      metal-nitride; Li5N3Ni3, ICSD #411152
A2BC_hP4_191_c_a_b-001.{ABC}                     hP4       P6/mmm              -         Li2N1Na1            #191      metal-nitride; Li2N1Na1, ICSD #92308
ABC2_hP4_191_b_a_c-001.{ABC}                     hP4       P6/mmm              -         Li1N1Na2            #191      metal-nitride; Li1N1Na2, ICSD #92310
ABC2_hP8_194_c_a_f-001.{ABC}                     hP8       P6_{3}/mmc          -         Ba1Ce1N2            #194      metal-nitride; Ba1Ce1N2, ICSD #74791
AB2C_hP8_194_a_f_b-001.{ABC}                     hP8       P6_{3}/mmc          -         Fe1N2W1             #194      metal-nitride; Fe1N2W1, ICSD #75971
ABC2_hP8_194_a_c_f-001.{ABC}                     hP8       P6_{3}/mmc          -         Mg1Mo1N2            #194      metal-nitride; Mg1Mo1N2, ICSD #185913
A3BC_hP10_194_h_a_c-001.{ABC}                    hP10      P6_{3}/mmc          -         Ba3N1Na1            #194      metal-nitride; Ba3N1Na1, ICSD #67497
AB3C4_hP16_194_c_ae_2f-001.{ABC}                 hP16      P6_{3}/mmc          -         Al1N3V4             #194      metal-nitride; Al1N3V4, ICSD #181351
AB2C_cF32_227_b_c_a.{ABC}                        cF32      Fd-3m               -         Cs1N2Nb1            #227      metal-nitride; Cs1N2Nb1, ICSD #72546

*****************************************************************************************************************************************************
ANRL - metals - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC2_mC16_8_b_2a_2ab-001.{ABC}                   mC16      Cm                  -         Bi1Pb1Pd2           #8        metallic; Bi1Pb1Pd2, ICSD #58830
A2B4C3_mP9_10_m_2n_am-001.{ABC}                  mP9       P2/m                -         Co2Gd4Mg3           #10       metallic; Co2Gd4Mg3, ICSD #417035
A2B5C4_mP11_10_m_a2m_2n-001.{ABC}                mP11      P2/m                -         Co2In5Zr4           #10       metallic; Co2In5Zr4, ICSD #55578
AB2C2_mP10_11_e_2e_2e-001.ABC/BCA/CAB            mP10      P2_{1}/m            -         Ba1Bi2Pd2           #11       metallic; Ba1Bi2Pd2, ICSD #416299
A2BC2_mP10_11_2e_e_2e-001.ABC/BCA/CAB            mP10      P2_{1}/m            -         Sn2Sr1Zn2           #11       metallic; Sn2Sr1Zn2, ICSD #424108
A2B4C_mP14_11_2e_4e_e-001.{ABC}                  mP14      P2_{1}/m            -         Au2In4Yb1           #11       metallic; Au2In4Yb1, ICSD #261042
A2BC2_mC10_12_i_a_i-001.{ABC}                    mC10      C2/m                -         Cu2Eu1Sn2           #12       metallic; Cu2Eu1Sn2, ICSD #182050
ABC_mC12_12_i_i_i-001.{ABC}                      mC12      C2/m                -         Al1Ce1Co1           #12       metallic; Al1Ce1Co1, ICSD #20632
A3B2C2_mC14_12_ai_i_i-001.{ABC}                  mC14      C2/m                -         In3Rh2Sr2           #12       metallic; In3Rh2Sr2, ICSD #410985
A6B2C_mC18_12_3i_i_a-001.{ABC}                   mC18      C2/m                -         Na6Sn2Zn1           #12       metallic; Na6Sn2Zn1, ICSD #260159
A2B5C4_mC22_12_i_a2i_2i-001.{ABC}                mC22      C2/m                -         Mn2Sn5Yb4           #12       metallic; Mn2Sn5Yb4, ICSD #419134
A3B4C8_mC30_12_ai_hi_2ij-001.{ABC}               mC30      C2/m                -         Ba3Li4Sn8           #12       metallic; Ba3Li4Sn8, ICSD #240016
A6B2C7_mC30_12_3i_i_a3i-001.{ABC}                mC30      C2/m                -         Ca6Cu2Sn7           #12       metallic; Ca6Cu2Sn7, ICSD #171243
AB2C2_mC20_15_e_f_f-001.{ABC}                    mC20      C2/c                -         In1Pd2Sr2           #15       metallic; In1Pd2Sr2, ICSD #391432
A2B4C_mC28_15_f_2f_e-001.{ABC}                   mC28      C2/c                -         Na2Sn4Sr1           #15       metallic; Na2Sn4Sr1, ICSD #261117
ABC2_oC16_20_b_a_c-001.{ABC}                     oC16      C222_{1}            -         Ag1Bi1K2            #20       metallic; Ag1Bi1K2, ICSD #1156
ABC_oP12_33_a_a_a-001.ABC/BCA/CAB                oP12      Pna2_{1}            -         Ce1Pd1Sn1           #33       metallic; Ce1Pd1Sn1, ICSD #106418
ABC_oC12_36_a_a_a-001.ABC/BCA/CAB                oC12      Cmc2_{1}            -         Ag1Ce1Sn1           #36       metallic; Ag1Ce1Sn1, ICSD #55819
ABC2_oC32_36_b_2a_2ab-001.{ABC}                  oC32      Cmc2_{1}            -         Bi1Pb1Pd2           #36       metallic; Bi1Pb1Pd2, ICSD #56278
ABC_oC18_38_ad_be_bd-001.{ABC}                   oC18      Amm2                -         In1Nd1Rh1           #38       metallic; In1Nd1Rh1, ICSD #59430
A10BC3_oP14_47_2q2rs_h_bt-001.{ABC}              oP14      Pmmm                -         Ga10Ni1Pr3          #47       metallic; Ga10Ni1Pr3, ICSD #20663
ABC2_oP8_51_f_e_be-001.{ABC}                     oP8       Pmma                -         In1La1Ni2           #51       metallic; In1La1Ni2, ICSD #106800
A2BC_oP8_51_ae_f_f-001.{ABC}                     oP8       Pmma                -         Co2Ga1La1           #51       metallic; Co2Ga1La1, ICSD #623100
A4BC_oP12_51_afj_e_e-001.{ABC}                   oP12      Pmma                -         Al4Co1La1           #51       metallic; Al4Co1La1, ICSD #9986
A3B2C_oP12_51_ak_ef_f-001.{ABC}                  oP12      Pmma                -         Au3Rb2Tl1           #51       metallic; Au3Rb2Tl1, ICSD #249924
AB4C_oP12_51_e_afj_e-001.{ABC}                   oP12      Pmma                -         Ca1In4Rh1           #51       metallic; Ca1In4Rh1, ICSD #410891
ABC_oP12_57_c_d_d-001.{ABC}                      oP12      Pbcm                -         Al1Ca1Pd1           #57       metallic; Al1Ca1Pd1, ICSD #370036
AB2C_oP16_57_d_2d_c-001.{ABC}                    oP16      Pbcm                -         Bi1K2Sn1            #57       metallic; Bi1K2Sn1, ICSD #107616
AB2C_oP8_59_a_e_b-001.ABC/BCA/CAB                oP8       Pmmn                -         Cu1Ni2Ti1           #59       metallic; Cu1Ni2Ti1, ICSD #628580
A3B3C_oP14_59_ae_ae_b-001.{ABC}                  oP14      Pmmn                -         Au3In3Sr1           #59       metallic; Au3In3Sr1, ICSD #245679
AB3C4_oP16_59_b_ae_ef-001.{ABC}                  oP16      Pmmn                -         Au1K3Sn4            #59       metallic; Au1K3Sn4, ICSD #107444
AB3C4_oP16_59_b_ae_ef-002.{ABC}                  oP16      Pmmn                -         Au1Cs3Pb4           #59       metallic; Au1Cs3Pb4, ICSD #107448
ABC_oP12_62_c_c_c-003.{ABC}                      oP12      Pnma                -         Ca1Mg1Sn1           #62       metallic; Ca1Mg1Sn1, ICSD #42757, equivalent to CuMnP (http://aflow.org/prototype-encyclopedia/ABC_oP12_62_c_c_c.html, part 3)
ABC_oP12_62_c_c_c-004.ABC/BCA/CAB                oP12      Pnma                -         Ni1Sn1Tb1           #62       metallic; Ni1Sn1Tb1, ICSD #54301
ABC_oP12_62_c_c_c-005.{ABC}                      oP12      Pnma                -         Co1Dy1Sn1           #62       metallic; Co1Dy1Sn1, ICSD #54415
ABC_oP12_62_c_c_c-006.{ABC}                      oP12      Pnma                -         Bi1Ca1Li1           #62       metallic; Bi1Ca1Li1, ICSD #58762
ABC_oP12_62_c_c_c-007.{ABC}                      oP12      Pnma                -         Ni1Sc1Sn1           #62       metallic; Ni1Sc1Sn1, ICSD #105338
ABC_oP12_62_c_c_c-008.ABC/BCA/CAB                oP12      Pnma                -         Co1Dy1Sn1           #62       metallic; Co1Dy1Sn1, ICSD #106459
ABC_oP12_62_c_c_c-009.ABC/BCA/CAB                oP12      Pnma                -         La1Ni1Sn1           #62       metallic; La1Ni1Sn1, ICSD #108571
ABC_oP12_62_c_c_c-010.{ABC}                      oP12      Pnma                -         Al1Ce1Pt1           #62       metallic; Al1Ce1Pt1, ICSD #150172
ABC_oP12_62_c_c_c-011.{ABC}                      oP12      Pnma                -         La1Ni1Sn1           #62       metallic; La1Ni1Sn1, ICSD #157921
ABC_oP12_62_c_c_c-012.{ABC}                      oP12      Pnma                -         Pt1Yb1Zn1           #62       metallic; Pt1Yb1Zn1, ICSD #159305
ABC_oP12_62_c_c_c-013.{ABC}                      oP12      Pnma                -         Er1Ga1Rh1           #62       metallic; Er1Ga1Rh1, ICSD #630573
ABC_oP12_62_c_c_c-014.{ABC}                      oP12      Pnma                -         Ga1Pd1Sc1           #62       metallic; Ga1Pd1Sc1, ICSD #635080
AB2C_oP16_62_c_d_c-001.{ABC}                     oP16      Pnma                -         Ca1Cd2Pd1           #62       metallic; Ca1Cd2Pd1, ICSD #425509
AB2C_oP16_62_c_d_c-002.{ABC}                     oP16      Pnma                -         Ga1Pd2Y1            #62       metallic; Ga1Pd2Y1, ICSD #602790
ABC_oC12_63_c_c_c-001.{ABC}                      oC12      Cmcm                -         Mg1Ni1Tb1           #63       metallic; Mg1Ni1Tb1, ICSD #166874
ABC_oC12_63_c_c_a-001.{ABC}                      oC12      Cmcm                -         Pt1Sc1Zn1           #63       metallic; Pt1Sc1Zn1, ICSD #424600
A2BC_oC16_63_f_c_c-001.{ABC}                     oC16      Cmcm                -         Al2Co1Y1            #63       metallic; Al2Co1Y1, ICSD #57645
A2BC_oC16_63_f_c_c-002.{ABC}                     oC16      Cmcm                E1_{1}    Al2Cu1Mg1           #63       metallic; Al2Cu1Mg1, ICSD #57693, equivalent to Al2CuMg (E1_{a}, http://aflow.org/prototype-encyclopedia/A2BC_oC16_63_f_c_c.html, part 3)
ABC2_oC16_63_a_c_g-001.{ABC}                     oC16      Cmcm                -         Au1Bi1Li2           #63       metallic; Au1Bi1Li2, ICSD #261785
ABC2_oC16_63_a_c_g-002.{ABC}                     oC16      Cmcm                -         Au1Bi1Li2           #63       metallic; Au1Bi1Li2, ICSD #261786
ABC2_oC16_63_a_c_g-003.{ABC}                     oC16      Cmcm                -         Au1Bi1Na2           #63       metallic; Au1Bi1Na2, ICSD #261788
ABC2_oC16_63_a_c_g-004.{ABC}                     oC16      Cmcm                -         Au1Bi1K2            #63       metallic; Au1Bi1K2, ICSD #380341
ABC2_oC16_63_c_c_2c-001.{ABC}                    oC16      Cmcm                -         Ba1Cu1Sn2           #63       metallic; Ba1Cu1Sn2, ICSD #58647
ABC2_oC16_63_c_c_2c-002.{ABC}                    oC16      Cmcm                -         Co1Ho1Sn2           #63       metallic; Co1Ho1Sn2, ICSD #240098
ABC2_oC16_63_c_c_2c-003.{ABC}                    oC16      Cmcm                -         Ce1Ni1Sn2           #63       metallic; Ce1Ni1Sn2, ICSD #621687
AB2C_oC16_63_c_2c_c-003.{ABC}                    oC16      Cmcm                -         Co1Sn2Tb1           #63       metallic; Co1Sn2Tb1, ICSD #240096
AB2C_oC16_63_c_f_c-001.{ABC}                     oC16      Cmcm                -         La1Mg2Ni1           #63       metallic; La1Mg2Ni1, ICSD #96152
A2BC_oC16_63_g_c_c-001.{ABC}                     oC16      Cmcm                -         Cd2Cu1Er1           #63       metallic; Cd2Cu1Er1, ICSD #99139
A2BC_oC16_63_f_c_c-003.{ABC}                     oC16      Cmcm                -         Ga2Ni1Sc1           #63       metallic; Ga2Ni1Sc1, ICSD #103875
A2BC_oC16_63_f_c_c-004.{ABC}                     oC16      Cmcm                -         Mg2Pd1Sr1           #63       metallic; Mg2Pd1Sr1, ICSD #425493
A2BC_oC16_63_f_c_c-005.{ABC}                     oC16      Cmcm                -         Al2Ce1Pt1           #63       metallic; Al2Ce1Pt1, ICSD #658140
A2B2C_oC20_63_g_2c_a-001.{ABC}                   oC20      Cmcm                -         Bi2Cs2Pt1           #63       metallic; Bi2Cs2Pt1, ICSD #658701
A4BC_oC24_63_acf_c_c-001.{ABC}                   oC24      Cmcm                -         Al4Dy1Ni1           #63       metallic; Al4Dy1Ni1, ICSD #57760
AB3C2_oC24_63_c_cg_e-001.{ABC}                   oC24      Cmcm                -         Ca1Ga3Ni2           #63       metallic; Ca1Ga3Ni2, ICSD #58898
A3BC3_oC28_63_cf_a_cf-001.{ABC}                  oC28      Cmcm                -         Co3Ga1Y3            #63       metallic; Co3Ga1Y3, ICSD #10044
ABC2_oC32_63_ac_f_2cf-001.{ABC}                  oC32      Cmcm                -         La1Rh1Sn2           #63       metallic; La1Rh1Sn2, ICSD #410732
A2B2C_oC10_65_j_i_a-001.{ABC}                    oC10      Cmmm                -         Ho2Ni2Pb1           #65       metallic; Ho2Ni2Pb1, ICSD #54612
A2B3C_oC12_65_i_cf_a-001.{ABC}                   oC12      Cmmm                -         Al2Ni3Pr1           #65       metallic; Al2Ni3Pr1, ICSD #107864
A3B2C2_oC14_65_aj_j_i-001.{ABC}                  oC14      Cmmm                -         Mg3Ni2Tb2           #65       metallic; Mg3Ni2Tb2, ICSD #240761
A2BC_oC16_65_aci_j_i-001.{ABC}                   oC16      Cmmm                -         Ga2Nd1Ni1           #65       metallic; Ga2Nd1Ni1, ICSD #103850
A2BC6_oC18_65_j_a_2ij-001.{ABC}                  oC18      Cmmm                -         Dy2Ni1Sn6           #65       metallic; Dy2Ni1Sn6, ICSD #165783
A3B2C7_oC24_65_aj_i_c2ij-001.{ABC}               oC24      Cmmm                -         Ce3Ni2Sn7           #65       metallic; Ce3Ni2Sn7, ICSD #102239
A3B2C7_oC24_65_aj_i_c2ij-002.{ABC}               oC24      Cmmm                -         La3Ni2Sn7           #65       metallic; La3Ni2Sn7, ICSD #160913
A2B2C_oI10_71_f_h_a-001.{ABC}                    oI10      Immm                -         Ce2Ni2Sn1           #71       metallic; Ce2Ni2Sn1, ICSD #55495
A2B2C_oI10_71_f_h_a-002.{ABC}                    oI10      Immm                -         Ca2Cu2Ga1           #71       metallic; Ca2Cu2Ga1, ICSD #58885
A2B2C3_oI14_71_e_f_af-001.{ABC}                  oI14      Immm                -         Al2Sn2Sr3           #71       metallic; Al2Sn2Sr3, ICSD #9564
A5B2C_oI16_71_an_e_d-001.{ABC}                   oI16      Immm                -         Al5Ni2Pr1           #71       metallic; Al5Ni2Pr1, ICSD #58046
A4B3C4_oI22_71_n_af_eg-001.{ABC}                 oI22      Immm                -         Ag8Ce6Sn8           #71       metallic; Ag8Ce6Sn8, ICSD #55826
A4B3C4_oI22_71_n_af_eg-002.{ABC}                 oI22      Immm                -         Ag8Pr6Sn8           #71       metallic; Ag8Pr6Sn8, ICSD #55827
A4B3C4_oI22_71_n_af_eh-001.{ABC}                 oI22      Immm                -         Ag4Dy3Sn4           #71       metallic; Ag4Dy3Sn4, ICSD #156968
A3B4C4_oI22_71_af_eh_n-001.{ABC}                 oI22      Immm                -         La3Pd4Zn4           #71       metallic; La3Pd4Zn4, ICSD #182774
A6B6C_oI26_71_hk_efg_a-001.{ABC}                 oI26      Immm                -         Fe6Ga6Sc1           #71       metallic; Fe6Ga6Sc1, ICSD #103465
A3B3C8_oI28_71_cg_be_2n-001.{ABC}                oI28      Immm                -         Ag3Er3Ga8           #71       metallic; Ag3Er3Ga8, ICSD #605107
A3B3C8_oI28_71_cg_be_2n-002.{ABC}                oI28      Immm                -         Au3Er3Ga8           #71       metallic; Au3Er3Ga8, ICSD #611824
A3B10C_oI28_71_be_g2n_c-001.{ABC}                oI28      Immm                -         Ba3Hg10In1          #71       metallic; Ba3Hg10In1, ICSD #290305
A6B5C3_oI28_71_gn_cn_be-001.{ABC}                oI28      Immm                -         Cu6Sn5Yb3           #71       metallic; Cu6Sn5Yb3, ICSD #413485
A9B3C2_oI28_71_c2n_be_g-001.{ABC}                oI28      Immm                -         Ga9Ho3Pd2           #71       metallic; Ga9Ho3Pd2, ICSD #634397
A2B2C_oI20_72_j_j_a-001.{ABC}                    oI20      Ibam                -         Ba2Bi2Zn1           #72       metallic; Ba2Bi2Zn1, ICSD #421424
ABC3_tI10_107_a_a_ab-001.{ABC}                   tI10      I4mm                -         Ba1Ni1Sn3           #107      metallic; Ba1Ni1Sn3, ICSD #58662, equivalent to BaNiSn3 (http://aflow.org/prototype-encyclopedia/ABC3_tI10_107_a_a_ab.html, part 3)
A3BC_tI10_107_ab_a_a-001.{ABC}                   tI10      I4mm                -         Al3Cu1Pr1           #107      metallic; Al3Cu1Pr1, ICSD #290390
A3BC_tI10_119_bf_a_c-001.{ABC}                   tI10      I-4m2               -         In3Mg1Sr1           #119      metallic; In3Mg1Sr1, ICSD #249592
AB4C2_tI28_120_d_i_e-001.{ABC}                   tI28      I-4c2               -         Ce1Ni4Sn2           #120      metallic; Ce1Ni4Sn2, ICSD #102238
A4B6C_tI22_121_i_ci_a-001.{ABC}                  tI22      I-42m               -         Ru4Sn6Y1            #121      metallic; Ru4Sn6Y1, ICSD #54354
ABC2_tP4_123_a_c_e-002.ABC/BCA/CAB               tP4       P4/mmm              -         Fe1Ni1Pt2           #123      metallic; Fe1Ni1Pt2, ICSD #42564
AB2C_tP4_123_b_h_a-001.ABC/BCA/CAB               tP4       P4/mmm              -         Cd1Pt2Zn1           #123      metallic; Cd1Pt2Zn1, ICSD #102057
AB4C_tP6_123_b_i_a-001.{ABC}                     tP6       P4/mmm              -         Co1Ga4Hf1           #123      metallic; Co1Ga4Hf1, ICSD #623081
AB5C_tP7_123_a_ci_b-001.{ABC}                    tP7       P4/mmm              -         Co1In5Tb1           #123      metallic; Co1In5Tb1, ICSD #623944
AB6C_tP8_123_a_hi_b-001.{ABC}                    tP8       P4/mmm              -         Ce1Ga6Pd1           #123      metallic; Ce1Ga6Pd1, ICSD #240161
AB8C2_tP11_123_a_ehi_g-001.{ABC}                 tP11      P4/mmm              -         Co1Ga8Ho2           #123      metallic; Co1Ga8Ho2, ICSD #42426
AB8C2_tP11_123_a_ehi_g-002.{ABC}                 tP11      P4/mmm              -         Fe1Ga8Ho2           #123      metallic; Fe1Ga8Ho2, ICSD #180132
A2BC2_tP10_127_g_a_h-001.{ABC}                   tP10      P4/mbm              -         Ni2Sn1Zr2           #127      metallic; Ni2Sn1Zr2, ICSD #54303
AB2C2_tP10_127_a_g_h-001.{ABC}                   tP10      P4/mbm              -         Al1Cu2Re2           #127      metallic; Al1Cu2Re2, ICSD #57706
A2BC_tP8_129_bc_c_a-001.{ABC}                    tP8       P4/nmm              -         Sn2Tb1Zn1           #129      metallic; Sn2Tb1Zn1, ICSD #163425
A2BC_tP8_129_bc_c_a-002.{ABC}                    tP8       P4/nmm              -         Bi2La1Li1           #129      metallic; Bi2La1Li1, ICSD #415728
AB2C2_tP10_129_c_ac_bc-001.ABC/BCA/CAB           tP10      P4/nmm              -         Nd1Ni2Sn2           #129      metallic; Nd1Ni2Sn2, ICSD #160053
AB2C2_tP10_129_c_bc_ac-001.ABC/BCA/CAB           tP10      P4/nmm              -         Eu1Sn2Zn2           #129      metallic; Eu1Sn2Zn2, ICSD #162267, equivalent to Be2CaGe2 (http://aflow.org/prototype-encyclopedia/A2BC2_tP10_129_ac_c_bc.html, part 3)
A2BC2_tP10_129_bc_c_ac-001.ABC/BCA/CAB           tP10      P4/nmm              -         Bi2Ce1Ni2           #129      metallic; Bi2Ce1Ni2, ICSD #616562
A2BC4_tP14_129_ac_b_j-001.{ABC}                  tP14      P4/nmm              -         Ce2Ru1Zn4           #129      metallic; Ce2Ru1Zn4, ICSD #418547
AB2C_tI8_139_b_d_a-001.ABC/BCA/CAB               tI8       I4/mmm              -         Li1Pd2Tl1           #139      metallic; Li1Pd2Tl1, ICSD #54365
A2BC_tI8_139_d_a_b-001.ABC/BCA/CAB               tI8       I4/mmm              -         Co2Ga1Ni1           #139      metallic; Co2Ga1Ni1, ICSD #169733
AB2C2_tI10_139_a_d_e-001.{ABC}                   tI10      I4/mmm              -         Ba1Mn2Sn2           #139      metallic; Ba1Mn2Sn2, ICSD #405
A2BC2_tI10_139_d_a_e-001.{ABC}                   tI10      I4/mmm              -         Ag2Ba1Sn2           #139      metallic; Ag2Ba1Sn2, ICSD #25332
A2BC2_tI10_139_d_a_e-002.{ABC}                   tI10      I4/mmm              -         Al2Ce1Ga2           #139      metallic; Al2Ce1Ga2, ICSD #55789
A2BC2_tI10_139_d_a_e-003.{ABC}                   tI10      I4/mmm              -         Al2Ca1Zn2           #139      metallic; Al2Ca1Zn2, ICSD #57550
ABC_tI12_139_c_e_e-001.{ABC}                     tI12      I4/mmm              -         Mg1Pr1Sn1           #139      metallic; Mg1Pr1Sn1, ICSD #182479
ABC_tI12_139_e_c_e-001.ABC/BCA/CAB               tI12      I4/mmm              -         Au1Cu1Sn1           #139      metallic; Au1Cu1Sn1, ICSD #611768
A4B2C_tI14_139_h_d_a-001.{ABC}                   tI14      I4/mmm              -         Al4Mo2Yb1           #139      metallic; Al4Mo2Yb1, ICSD #456
AB4C2_tI14_139_a_h_d-001.{ABC}                   tI14      I4/mmm              -         Er1Ga4Ti2           #139      metallic; Er1Ga4Ti2, ICSD #630579
A5B2C_tI16_139_bg_e_a-001.{ABC}                  tI16      I4/mmm              -         Al5Ni2Zr1           #139      metallic; Al5Ni2Zr1, ICSD #58084
A2BC_tI16_139_ce_e_d-001.{ABC}                   tI16      I4/mmm              -         Bi2Sr1Zn1           #139      metallic; Bi2Sr1Zn1, ICSD #41924
AB2C_tI16_139_e_ce_d-001.{ABC}                   tI16      I4/mmm              -         Ba1Bi2Cd1           #139      metallic; Ba1Bi2Cd1, ICSD #58635
AB2C_tI16_139_e_ce_d-002.{ABC}                   tI16      I4/mmm              -         Ba1Bi2Zn1           #139      metallic; Ba1Bi2Zn1, ICSD #58638
A2B5C3_tI20_139_e_bg_ae-001.{ABC}                tI20      I4/mmm              -         K2Mg5Sn3            #139      metallic; K2Mg5Sn3, ICSD #421342
A8BC4_tI26_139_ij_a_f-001.{ABC}                  tI26      I4/mmm              -         Al8Ca1Cu4           #139      metallic; Al8Ca1Cu4, ICSD #57539
A4B8C_tI26_139_i_fj_a-001.{ABC}                  tI26      I4/mmm              -         Cr4Fe8Y1            #139      metallic; Cr4Fe8Y1, ICSD #168240
A8B4C_tI26_139_ij_f_a-001.{ABC}                  tI26      I4/mmm              -         Al8Mn4Y1            #139      metallic; Al8Mn4Y1, ICSD #57997
ABC4_tI24_140_a_c_l-001.{ABC}                    tI24      I4/mcm              -         Ir1Li1Sn4           #140      metallic; Ir1Li1Sn4, ICSD #172149
A2BC4_tI28_140_h_a_k-001.{ABC}                   tI28      I4/mcm              -         Bi2Mn1Ti4           #140      metallic; Bi2Mn1Ti4, ICSD #150145
A4B2C_tI28_140_l_h_a-001.{ABC}                   tI28      I4/mcm              -         Cu4Sn2Sr1           #140      metallic; Cu4Sn2Sr1, ICSD #182104, equivalent to Sb2SiV4 (http://aflow.org/prototype-encyclopedia/A2BC4_tI28_140_h_a_k.html, part 3)
AB5C2_tI32_140_a_cl_h-001.{ABC}                  tI32      I4/mcm              -         Bi1Er5Pt2           #140      metallic; Bi1Er5Pt2, ICSD #107217
A3B4C_tI32_140_ah_l_c-001.{ABC}                  tI32      I4/mcm              -         Bi3In4Pb1           #140      metallic; Bi3In4Pb1, ICSD #616725
AB5C2_tI32_140_a_bk_h-001.{ABC}                  tI32      I4/mcm              -         Ga1Nb5Sn2           #140      metallic; Ga1Nb5Sn2, ICSD #103842
A2BC5_tI32_140_h_a_cl-001.{ABC}                  tI32      I4/mcm              -         Au2Bi1Tb5           #140      metallic; Au2Bi1Tb5, ICSD #156956
AB2C5_tI32_140_a_h_cl-001.{ABC}                  tI32      I4/mcm              -         Bi1Co2Ho5           #140      metallic; Bi1Co2Ho5, ICSD #161658
A4B7C3_hR14_146_ab_a2b_b-001.{ABC}               hR14      R3                  -         Ga4Pd7Zn3           #146      metallic; Ga4Pd7Zn3, ICSD #103911
ABC_hP9_156_3a_b2c_2bc-001.{ABC}                 hP9       P3m1                -         Ca1Li1Sn1           #156      metallic; Ca1Li1Sn1, ICSD #58911
A5BC4_hR10_160_5a_a_4a-001.{ABC}                 hR10      R3m                 -         Li5Na1Sn4           #160      metallic; Li5Na1Sn4, ICSD #12142
A2BC3_hP6_164_d_a_bd-002.{ABC}                   hP6       P-3m1               -         Bi2La1Li3           #164      metallic; Bi2La1Li3, ICSD #616769
A2BC_hR4_166_c_b_a-001.{ABC}                     hR4       R-3m                -         Cu2Ga1Sr1           #166      metallic; Cu2Ga1Sr1, ICSD #102938
AB2C_hR4_166_a_c_b-002.{ABC}                     hR4       R-3m                -         Ba1Cu2Ga1           #166      metallic; Ba1Cu2Ga1, ICSD #615828
AB2C4_hR7_166_a_c_2c-001.{ABC}                   hR7       R-3m                -         Au1Ni2Sn4           #166      metallic; Au1Ni2Sn4, ICSD #150127
AB4C6_hR11_166_a_2c_h-001.{ABC}                  hR11      R-3m                -         Gd1Sn4Ti6           #166      metallic; Gd1Sn4Ti6, ICSD #166193
A2BC_hR12_166_h_bc_ac-001.{ABC}                  hR12      R-3m                -         Al2Cu1Yb1           #166      metallic; Al2Cu1Yb1, ICSD #604213
A7B3C2_hR12_166_bh_ac_c-001.{ABC}                hR12      R-3m                -         Al7Ca3Cu2           #166      metallic; Al7Ca3Cu2, ICSD #57538
AB2C9_hR12_166_a_c_bch-001.{ABC}                 hR12      R-3m                -         La1Mg2Ni9           #166      metallic; La1Mg2Ni9, ICSD #55614
A9B2C_hR12_166_eh_c_a-001.{ABC}                  hR12      R-3m                -         Ni3Ti0.67Zr0.33     #166      metallic; Ni3Ti0.67Zr0.33, ICSD #646978
AB2C_hR12_166_bc_h_ac-001.{ABC}                  hR12      R-3m                -         Ag1Al2Pr1           #166      metallic; Ag1Al2Pr1, ICSD #604688
A2B7C3_hR12_166_c_bh_ac-001.{ABC}                hR12      R-3m                -         Ag2Al7Ca3           #166      metallic; Ag2Al7Ca3, ICSD #104173
A7B4C2_hR13_166_ah_2c_c-001.{ABC}                hR13      R-3m                -         Au7Rb4Sn2           #166      metallic; Au7Rb4Sn2, ICSD #58581
A3B2C2_hR14_166_3c_2c_abc-001.{ABC}              hR14      R-3m                -         Au3Sn2Yb2           #166      metallic; Au3Sn2Yb2, ICSD #710044
ABC_hP6_186_b_b_a-002.ABC/BCA/CAB                hP6       P6_{3}mc            -         Cu1Sn1Ti1           #186      metallic; Cu1Sn1Ti1, ICSD #54657
ABC_hP6_186_b_a_b-001.ABC/BCA/CAB                hP6       P6_{3}mc            -         Au1Pr1Sn1           #186      metallic; Au1Pr1Sn1, ICSD #54997
ABC_hP6_186_b_a_b-002.{ABC}                      hP6       P6_{3}mc            -         Bi1Li1Zn1           #186      metallic; Bi1Li1Zn1, ICSD #100115
ABC_hP6_186_a_b_b-001.ABC/BCA/CAB                hP6       P6_{3}mc            -         Ce1Cu1Sn1           #186      metallic; Ce1Cu1Sn1, ICSD #156394
A4BC_hP12_186_ac_b_b-001.ABC/BCA/CAB             hP12      P6_{3}mc            -         Cu4In1Mn1           #186      metallic; Cu4In1Mn1, ICSD #424277
A3BC4_hP16_186_a2b_b_2a2b-001.{ABC}              hP16      P6_{3}mc            -         Au3Li1Sn4           #186      metallic; Au3Li1Sn4, ICSD #412207
ABC_hP6_187_i_h_ab-001.ABC/BCA/CAB               hP6       P-6m2               -         Au1Ga1Zr1           #187      metallic; Au1Ga1Zr1, ICSD #156264
A4BC_hP6_187_hi_b_a-001.ABC/BCA/CAB              hP6       P-6m2               -         Ga4Li1Y1            #187      metallic; Ga4Li1Y1, ICSD #98666
ABC_hP9_187_dh_ai_eg-001.{ABC}                   hP9       P-6m2               -         Ca1Li1Pb1           #187      metallic; Ca1Li1Pb1, ICSD #409533
A5B5C_hP11_187_hk_cgi_a-001.{ABC}                hP11      P-6m2               -         Al5Ba5Pb1           #187      metallic; Al5Ba5Pb1, ICSD #416337
A4B6C_hP11_187_ck_jk_a-001.{ABC}                 hP11      P-6m2               -         Au4In6K1            #187      metallic; Au4In6K1, ICSD #249520
A3B2C3_hP8_189_f_c_g-001.{ABC}                   hP8       P-62m               -         In3Rh2Ti3           #189      metallic; In3Rh2Ti3, ICSD #410967
ABC_hP9_189_g_f_ad-001.{ABC}                     hP9       P-62m               -         Bi1Dy1Rh1           #189      metallic; Bi1Dy1Rh1, ICSD #51845, equivalent to ZlNiZr (http://aflow.org/prototype-encyclopedia/ABC_hP9_189_g_ad_f.html, part 3)
AB2C6_hP9_189_b_c_fg-001.{ABC}                   hP9       P-62m               -         Co1Ga2Zr6           #189      metallic; Co1Ga2Zr6, ICSD #20876
A5B3C_hP9_189_cf_g_b-001.{ABC}                   hP9       P-62m               -         Ga5Lu3Ni1           #189      metallic; Ga5Lu3Ni1, ICSD #634575
ABC_hP9_189_g_bc_f-001.{ABC}                     hP9       P-62m               -         Pb1Pd1Y1            #189      metallic; Pb1Pd1Y1, ICSD #54315
A2BC6_hP9_189_c_b_fg-001.{ABC}                   hP9       P-62m               -         Bi2Fe1Ho6           #189      metallic; Bi2Fe1Ho6, ICSD #96253
ABC_hP9_189_f_bc_g-001.{ABC}                     hP9       P-62m               -         Mg1Sn1Yb1           #189      metallic; Mg1Sn1Yb1, ICSD #54344
ABC_hP9_189_f_bc_g-002.{ABC}                     hP9       P-62m               -         Ga1Ni1Ti1           #189      metallic; Ga1Ni1Ti1, ICSD #103885
ABC_hP9_189_f_g_bc-001.{ABC}                     hP9       P-62m               -         Al1Dy1Ni1           #189      metallic; Al1Dy1Ni1, ICSD #107416
ABC_hP9_189_f_g_bc-002.{ABC}                     hP9       P-62m               -         Al1Hf1Pt1           #189      metallic; Al1Hf1Pt1, ICSD #608141
ABC_hP9_189_f_g_bc-003.{ABC}                     hP9       P-62m               -         Al1Nd1Ni1           #189      metallic; Al1Nd1Ni1, ICSD #608770
ABC_hP9_189_bc_f_g-001.{ABC}                     hP9       P-62m               -         Cu1In1Tb1           #189      metallic; Cu1In1Tb1, ICSD #628137
A3B2C_hP6_191_g_c_a-001.{ABC}                    hP6       P6/mmm              -         Ag3Al2La1           #191      metallic; Ag3Al2La1, ICSD #57329
A3B2C_hP6_191_g_c_a-002.{ABC}                    hP6       P6/mmm              -         Co3Ga2Tb1           #191      metallic; Co3Ga2Tb1, ICSD #102452
A2BC3_hP6_191_c_a_g-001.{ABC}                    hP6       P6/mmm              -         Cu2Dy1In3           #191      metallic; Cu2Dy1In3, ICSD #108369
A2BC3_hP6_191_c_a_g-002.{ABC}                    hP6       P6/mmm              -         Al2Ce1Pt3           #191      metallic; Al2Ce1Pt3, ICSD #658142
A9BC2_hP12_191_fm_a_c-001.{ABC}                  hP12      P6/mmm              -         Al9Ba1Fe2           #191      metallic; Al9Ba1Fe2, ICSD #57518
A6BC6_hP13_191_i_a_cde-001.{ABC}                 hP13      P6/mmm              -         Mn6Sc1Sn6           #191      metallic; Mn6Sc1Sn6, ICSD #54273
A6B6C_hP13_191_i_cde_a-001.{ABC}                 hP13      P6/mmm              -         Mn6Sn6Tb1           #191      metallic; Mn6Sn6Tb1, ICSD #54277
A3B3C2_hP16_193_g_g_d-001.{ABC}                  hP16      P6_{3}/mcm          -         Ga3Hf3Nb2           #193      metallic; Ga3Hf3Nb2, ICSD #103733
ABC_hP6_194_d_a_c-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Pb1Sr1Zn1           #194      metallic; Pb1Sr1Zn1, ICSD #54319
ABC_hP6_194_c_a_d-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Bi1Ca1Cu1           #194      metallic; Bi1Ca1Cu1, ICSD #57018
ABC_hP6_194_c_a_d-002.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Cu1Gd1Sn1           #194      metallic; Cu1Gd1Sn1, ICSD #600994
ABC_hP6_194_c_d_a-002.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Al1Au1Ti1           #194      metallic; Al1Au1Ti1, ICSD #57506
ABC_hP6_194_a_c_d-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ba1Pb1Zn1           #194      metallic; Ba1Pb1Zn1, ICSD #106315
ABC_hP6_194_a_d_c-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ca1Hg1Pb1           #194      metallic; Ca1Hg1Pb1, ICSD #106345
AB2C_hP8_194_b_f_c-001.{ABC}                     hP8       P6_{3}/mmc          -         Al1Pt2Zr1           #194      metallic; Al1Pt2Zr1, ICSD #58139
AB2C_hP8_194_a_f_c-001.{ABC}                     hP8       P6_{3}/mmc          -         In1Pt2Y1            #194      metallic; In1Pt2Y1, ICSD #59505
A2BC_hP8_194_f_a_c-001.{ABC}                     hP8       P6_{3}/mmc          -         Cu2Li1Sn1           #194      metallic; Cu2Li1Sn1, ICSD #150602
A3BC2_hP12_194_h_a_f-001.{ABC}                   hP12      P6_{3}/mmc          -         Al3Ru1Sc2           #194      metallic; Al3Ru1Sc2, ICSD #58159
ABC_hP12_194_f_f_ab-001.{ABC}                    hP12      P6_{3}/mmc          -         Ga1Sn1Sr1           #194      metallic; Ga1Sn1Sr1, ICSD #66003
ABC_hP12_194_ab_f_f-001.{ABC}                    hP12      P6_{3}/mmc          -         Gd1Sn1Zn1           #194      metallic; Gd1Sn1Zn1, ICSD #152621
ABC_hP12_194_ab_f_f-002.{ABC}                    hP12      P6_{3}/mmc          -         La1Sn1Zn1           #194      metallic; La1Sn1Zn1, ICSD #152623
A3B3C2_hP16_194_h_af_f-001.{ABC}                 hP16      P6_{3}/mmc          -         Al2Ba2Ga1.33        #194      metallic; Al2Ba2Ga1.33, ICSD #415931
AB2C_cP16_215_e_bce_ad-001.{ABC}                 cP16      P-43m               -         Ir1Li2Mg1           #215      metallic; Ir1Li2Mg1, ICSD #180119
ABC2_cF16_216_a_d_bc.ABC/CAB/BCA                 cF16      F-43m               -         Ag1Al1Li2           #216      metallic; Ag1Al1Li2, ICSD #57330, equivalent to inverse Heusler CuHg2Ti (http://aflow.org/prototype-encyclopedia/AB2C_cF16_216_b_ad_c.html, part 3)
AB4C_cF24_216_a_e_c-001.ABC/BCA/CAB              cF24      F-43m               -         In1Ni4Zr1           #216      metallic; In1Ni4Zr1, ICSD #59462
A3B4C4_cI22_217_b_c_c-001.{ABC}                  cI22      I-43m               -         Ca3Ga4Ni4           #217      metallic; Ca3Ga4Ni4, ICSD #58899
A3BC_cP15_221_ag_c_d-001.{ABC}                   cP15      Pm-3m               -         Ga3Os1Tb1           #221      metallic; Ga3Os1Tb1, ICSD #412144

*****************************************************************************************************************************************************
ANRL - metal borides - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A4B_mC10_12_j_a-001.{AB}                         mC10      C2/m                -         B4Mn1               #12       metal-boride; B4Mn1, ICSD #15079 equivalent to 1034.AB,1034.BA [HTQC]
A2B_mC12_12_2i_i-002.{AB}                        mC12      C2/m                -         B2Mo1               #12       metal-boride; B2Mo1, ICSD #418398 equivalent to 1041.AB,1041.BA [HTQC]
A5B2_mC14_12_a2i_i-002.{AB}                      mC14      C2/m                -         B5W2                #12       metal-boride; B5W2, ICSD #20326 equivalent to 1019.AB,1019.BA,1046.AB,1046.BA [HTQC]
A3B4_mC28_15_ef_2f-001.{AB}                      mC28      C2/c                -         B3Ni4               #15       metal-boride; B3Ni4, ICSD #24308 equivalent to 1071.AB,1071.BA [HTQC]
A3B4_mC28_15_ef_2f-002.{AB}                      mC28      C2/c                -         B3Ni4               #15       metal-boride; B3Ni4, ICSD #150561 equivalent to 1070.AB,1070.BA [HTQC]
A4B_oP10_58_2g_a-001.{AB}                        oP10      Pnnm                -         B4Cr1               #58       metal-boride; B4Cr1, ICSD #186851 equivalent to 1055.AB,1055.BA,1056.AB,1056.BA [HTQC]
A2B_oP6_59_e_a-001.{AB}                          oP6       Pmmn                -         B2Ru1               #59       metal-boride; B2Ru1, ICSD #31871 equivalent to 1018.AB,1018.BA [HTQC], equivalent to B2Ru1 (http://aflow.org/prototype-encyclopedia/A2B_oP6_59_f_a.html, part 3)
AB_oP8_62_c_c-007.{AB}                           oP8       Pnma                -         B1Ti1               #62       metal-boride; B1Ti1, ICSD #24701 equivalent to 256,1020.AB,1020.BA,1050.AB,1050.BA [HTQC]
AB_oP8_62_c_c-008.{AB}                           oP8       Pnma                -         B1Li1               #62       metal-boride; B1Li1, ICSD #153291 equivalent to 1049.AB,1049.BA [HTQC]
AB_oP8_62_c_c-009.{AB}                           oP8       Pnma                -         B1Co1               #62       metal-boride; B1Co1, ICSD #612863 equivalent to 256,1020.AB,1020.BA,1047.AB,1047.BA,1050.AB,1050.BA [HTQC]
A2B_oP12_62_d_c-001.{AB}                         oP12      Pnma                -         B2Fe1               #62       metal-boride; B2Fe1, ICSD #425310 equivalent to 1066.AB,1066.BA [HTQC]
AB_oC8_63_a_c-001.{AB}                           oC8       Cmcm                -         B1Rh1               #63       metal-boride; B1Rh1, ICSD #150732 equivalent to 1029.AB,1029.BA [HTQC]
AB3_oC16_63_c_cf-001.{AB}                        oC16      Cmcm                -         B1Re3               #63       metal-boride; B1Re3, ICSD #43662 equivalent to 1022.AB,1022.BA [HTQC], equivalent to BRe3 (http://aflow.org/prototype-encyclopedia/AB3_oC16_63_c_cf.html, part 3)
A3B2_oC20_63_3c_2c-001.{AB}                      oC20      Cmcm                -         B3V2                #63       metal-boride; B3V2, ICSD #79258 equivalent to 1012.AB,1012.BA [HTQC]
A6B5_oC22_65_i2j_aij-001.{AB}                    oC22      Cmmm                -         B6Ta5               #65       metal-boride; B6Ta5, ICSD #68538 equivalent to 1058.AB,1058.BA [HTQC]
A4B_oI10_71_n_a-001.{AB}                         oI10      Immm                -         B4Cr1               #71       metal-boride; B4Cr1, ICSD #24353 equivalent to 1033.AB,1033.BA [HTQC]
A4B3_oI14_71_ef_af-001.{AB}                      oI14      Immm                -         B4Nb3               #71       metal-boride; B4Nb3, ICSD #76631 equivalent to 1045.AB,1045.BA [HTQC]
A4B3_oI14_71_ef_af-002.{AB}                      oI14      Immm                -         B4Mn3               #71       metal-boride; B4Mn3, ICSD #614732
A11B_tI24_119_aegi_d-001.{AB}                    tI24      I-4m2               -         B11Li1              #119      metal-boride; B11Li1, ICSD #164844 equivalent to 1062.AB,1062.BA [HTQC]
A11B_tI24_119_aegi_d-002.{AB}                    tI24      I-4m2               -         B11Li1              #119      metal-boride; B11Li1, ICSD #164845 equivalent to 1063.AB,1063.BA [HTQC]
AB2_tI12_121_ab_i-001.{AB}                       tI12      I-42m               C17       B1Fe2               #121      metal-boride; B1Fe2, ICSD #16809 equivalent to 1037.AB,1037.BA [HTQC], equivalent to Fe2B (C17, obsolete, http://aflow.org/prototype-encyclopedia/AB2_tI12_121_ab_i.html)
A12B_tI26_139_jm_a-001.{AB}                      tI26      I4/mmm              -         B12Sc1              #139      metal-boride; B12Sc1, ICSD #615424 equivalent to 1068.AB,1068.BA [HTQC]
AB2_tI12_140_b_h-001.{AB}                        tI12      I4/mcm              -         B1Fe2               #140      metal-boride; B1Fe2, ICSD #160789 equivalent to 1036.AB,1036.BA [HTQC]
AB2_tI12_140_b_h-002.{AB}                        tI12      I4/mcm              -         B1Fe2               #140      metal-boride; B1Fe2, ICSD #160790
AB2_tI12_140_c_h-001.{AB}                        tI12      I4/mcm              -         B1Ti2               #140      metal-boride; B1Ti2, ICSD #189385 equivalent to 1040.AB,1040.BA [HTQC]
A4B5_hR9_160_ab_2ab-001.{AB}                     hR9       R3m                 -         B4Li5               #160      metal-boride; B4Li5, ICSD #23041
A11B_hR12_160_2a3b_a-001.{AB}                    hR12      R3m                 -         B11Li1              #160      metal-boride; B11Li1, ICSD #164842 equivalent to 1060.AB,1060.BA [HTQC]
AB_hR4_166_c_c-003.{AB}                          hR4       R-3m                -         B2Li2               #166      metal-boride; B2Li2, ICSD #1 equivalent to 537,538,1032.AB,1032.BA [HTQC]
A2B_hR6_166_2c_c-001.{AB}                        hR6       R-3m                -         B4Mo2               #166      metal-boride; B4Mo2, ICSD #39554 equivalent to 1043.AB,1043.BA,1078.AB,1078.BA [HTQC]
A5B2_hR7_166_a2c_c-002.{AB}                      hR7       R-3m                -         B5W2                #166      metal-boride; B5W2, ICSD #20326
A3B_hR8_166_f_c-001.{AB}                         hR8       R-3m                -         B3Mo1               #166      metal-boride; B3Mo1, ICSD #167734 equivalent to 1052.AB,1052.BA [HTQC]
AB3_hP8_182_c_g-002.{AB}                         hP8       P6_{3}22            -         B1Fe3               #182      metal-boride; B1Fe3, ICSD #184958 equivalent to 1048.AB,1048.BA,1262.AB,1262.BA,1263.AB,1263.BA [HTQC]
A2B_hP3_191_d_a-001.{AB}                         hP3       P6/mmm              -         B2Zr1               #191      metal-boride; B2Zr1, ICSD #169458 equivalent to 1011.AB,1011.BA [HTQC]
A2B_hP3_191_d_a-002.{AB}                         hP3       P6/mmm              -         B2Ca1               #191      metal-boride; B2Ca1, ICSD #186764
A2B_hP3_191_d_a-003.{AB}                         hP3       P6/mmm              -         B2Nb1               #191      metal-boride; B2Nb1, ICSD #614887 equivalent to 519,520,1028.AB,1028.BA [HTQC]
AB_hP4_194_a_c-002.{AB}                          hP4       P6_{3}/mmc          -         B1Pt1               #194      metal-boride; B1Pt1, ICSD #24363 equivalent to 1016.AB,1016.BA [HTQC]
A2B_hP6_194_f_c-001.{AB}                         hP6       P6_{3}/mmc          -         B2Re1               #194      metal-boride; B2Re1, ICSD #23871 equivalent to 694.AB,694.BA,1017.AB,1017.BA [HTQC]
AB_hP8_194_bc_f-001.{AB}                         hP8       P6_{3}/mmc          -         B2Li2               #194      metal-boride; B2Li2, ICSD #2 equivalent to 1051.AB,1051.BA [HTQC]
A3B_hP8_194_af_c-001.{AB}                        hP8       P6_{3}/mmc          -         B3Re1               #194      metal-boride; B3Re1, ICSD #24361 equivalent to 1053.AB,1053.BA [HTQC], equivalent to B3Re (http://aflow.org/prototype-encyclopedia/A3B_hP8_194_af_c.html, part 3) 
A3B2_hP10_194_cf_f-002.{AB}                      hP10      P6_{3}/mmc          -         B3Ru2               #194      metal-boride; B3Ru2, ICSD #23715 equivalent to 1023.AB,1023.BA [HTQC]
A4B_hP10_194_2f_c-001.{AB}                       hP10      P6_{3}/mmc          -         B4Mo1               #194      metal-boride; B4Mo1, ICSD #167735 equivalent to 1057.AB,1057.BA [HTQC]
A3B2_hP10_194_bf_f-001.{AB}                      hP10      P6_{3}/mmc          -         B3Ru2               #194      metal-boride; B3Ru2, ICSD #108082
A2B_hP12_194_bcf_f-001.{AB}                      hP12      P6_{3}/mmc          -         B2W1                #194      metal-boride; B2W1, ICSD #23716 equivalent to 1067.AB,1067.BA [HTQC]
A5B2_hP14_194_fg_e-001.{AB}                      hP14      P6_{3}/mmc          -         B5Mo2               #194      metal-boride; B5Mo2, ICSD #167733 equivalent to 1072.AB,1072.BA [HTQC]
A11B_cF48_216_aeg_b-001.{AB}                     cF48      F-43m               -         B11Li1              #216      metal-boride; B11Li1, ICSD #164841 equivalent to 1059.AB,1059.BA [HTQC]
A11B_cF48_216_aef_d-001.{AB}                     cF48      F-43m               -         B11Li1              #216      metal-boride; B11Li1, ICSD #164843 equivalent to 1061.AB,1061.BA [HTQC]

*****************************************************************************************************************************************************
ANRL - metal borides - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A4B2C_aP7_2_2i_i_a-001.{ABC}                     aP7       P-1                 -         B4Co2Mo1            #2        metal-boride; B4Co2Mo1, ICSD #44175
A5B2C3_mC20_5_b2c_c_ac-001.{ABC}                 mC20      C2                  -         B5Ca2Os3            #5        metal-boride; B5Ca2Os3, ICSD #59227
A3BC_mP10_11_3e_e_e-001.{ABC}                    mP10      P2_{1}/m            -         B3Er1Mo1            #11       metal-boride; B3Er1Mo1, ICSD #65932
A2BC3_mC12_12_g_a_df-001.{ABC}                   mC12      C2/m                -         B2Er1Ir3            #12       metal-boride; B2Er1Ir3, ICSD #44230
A4B2C5_mC22_12_2i_i_aj-001.{ABC}                 mC22      C2/m                -         B4La2Ni5            #12       metal-boride; B4La2Ni5, ICSD #63501
A4B2C5_mC22_12_2i_i_aj-002.{ABC}                 mC22      C2/m                -         B4La2Ni5            #12       metal-boride; B4La2Ni5, ICSD #170618
A2B2C_mC20_15_f_f_e-001.{ABC}                    mC20      C2/c                -         B2Ni2Tb1            #15       metal-boride; B2Ni2Tb1, ICSD #57003
ABC_oF48_43_b_b_b-001.{ABC}                      oF48      Fdd2                -         B1Cu1Ir1            #43       metal-boride; B1Cu1Ir1, ICSD #75029
AB4C_oP6_47_d_rt_a-001.ABC/BCA/CAB               oP6       Pmmm                -         Al1B4Yb1            #47       metal-boride; Al1B4Yb1, ICSD #181368
AB4C3_oP8_47_a_qr_dt-001.{ABC}                   oP8       Pmmm                -         Al1B4Cr3            #47       metal-boride; Al1B4Cr3, ICSD #20082
A3B4C_oP16_47_qrs_qrst_ad-001.{ABC}              oP16      Pmmm                -         B3Ir4Zn1            #47       metal-boride; B3Ir4Zn1, ICSD #71640
A2B2C3_oP14_55_h_h_ag-001.{ABC}                  oP14      Pbam                -         B2Li2Rh3            #55       metal-boride; B2Li2Rh3, ICSD #417442
AB2C_oP16_62_c_d_c-003.{ABC}                     oP16      Pnma                -         B1Co2Fe1            #62       metal-boride; B1Co2Fe1, ICSD #603579
A2BC_oP16_62_2c_c_c-001.{ABC}                    oP16      Pnma                -         B2Co1Nb1            #62       metal-boride; B2Co1Nb1, ICSD #41893
A2BC_oP16_62_d_c_c-001.{ABC}                     oP16      Pnma                -         B2Ce1Ru1            #62       metal-boride; B2Ce1Ru1, ICSD #612853
ABC_oC12_63_c_c_c-002.{ABC}                      oC12      Cmcm                -         Al1B1Mo1            #63       metal-boride; Al1B1Mo1, ICSD #16777
ABC_oC12_63_c_c_c-003.{ABC}                      oC12      Cmcm                -         B1Nb1Ni1            #63       metal-boride; B1Nb1Ni1, ICSD #43005
A3BC_oC20_63_3c_c_c-001.{ABC}                    oC20      Cmcm                -         B3Co1V1             #63       metal-boride; B3Co1V1, ICSD #44188
A3BC3_oC28_63_cf_a_cf-002.{ABC}                  oC28      Cmcm                -         B3Co1W3             #63       metal-boride; B3Co1W3, ICSD #25753
A4BC3_oC32_63_2cf_a_cf-001.{ABC}                 oC32      Cmcm                -         B4Mg1Os3            #63       metal-boride; B4Mg1Os3, ICSD #91243
AB2C2_oC10_65_a_i_j-001.{ABC}                    oC10      Cmmm                -         Al1B2Cr2            #65       metal-boride; Al1B2Cr2, ICSD #20083
A2BC3_oC12_65_i_a_cf-001.{ABC}                   oC12      Cmmm                -         B2Ce1Rh3            #65       metal-boride; B2Ce1Rh3, ICSD #40777
A2B2C3_oC14_65_j_j_ai-001.{ABC}                  oC14      Cmmm                -         Al2B2Ru3            #65       metal-boride; Al2B2Ru3, ICSD #43843
A6B2C3_oC22_65_gp_j_ah-001.{ABC}                 oC22      Cmmm                -         B6Er2Ni3            #65       metal-boride; B6Er2Ni3, ICSD #65737
ABC4_oC24_65_i_bc_jkm-001.{ABC}                  oC24      Cmmm                -         B1Nd1Ni4            #65       metal-boride; B1Nd1Ni4, ICSD #20882
AB4C_oC24_65_h_gip_j-001.{ABC}                   oC24      Cmmm                -         Al1B4Yb1            #65       metal-boride; Al1B4Yb1, ICSD #167533
A4B2C5_oF44_69_o_g_acf-001.{ABC}                 oF44      Fmmm                -         B4Ca2Rh5            #69       metal-boride; B4Ca2Rh5, ICSD #66768
A2B2C_oF40_70_g_e_a-001.{ABC}                    oF40      Fddd                -         B2Rh2Sr1            #70       metal-boride; B2Rh2Sr1, ICSD #8153
A2BC2_oF40_70_f_a_e-001.{ABC}                    oF40      Fddd                -         B2Nd1Os2            #70       metal-boride; B2Nd1Os2, ICSD #601369
ABC_oF48_70_e_f_c-001.{ABC}                      oF48      Fddd                -         B1Ir1Li1            #70       metal-boride; B1Ir1Li1, ICSD #75027
A2BC2_oI10_71_h_a_f-001.{ABC}                    oI10      Immm                -         B2Co1W2             #71       metal-boride; B2Co1W2, ICSD #16776
A4B2C_oI14_71_ef_f_a-001.{ABC}                   oI14      Immm                -         B4Fe2Mo1            #71       metal-boride; B4Fe2Mo1, ICSD #44293
A4B2C_oI14_71_ef_f_a-002.{ABC}                   oI14      Immm                -         B4Mn2W1             #71       metal-boride; B4Mn2W1, ICSD #614782
A6BC2_oI18_71_hn_a_f-001.{ABC}                   oI18      Immm                -         B6Ce1Cr2            #71       metal-boride; B6Ce1Cr2, ICSD #16203
A6BC2_oI18_71_hn_a_f-002.{ABC}                   oI18      Immm                -         B6Ce1Cr2            #71       metal-boride; B6Ce1Cr2, ICSD #81544
ABC3_tP5_99_a_b_ac-001.{ABC}                     tP5       P4mm                -         B1Ce1Pt3            #99       metal-boride; B1Ce1Pt3, ICSD #95049
A2BC2_tI10_139_e_a_d-001.{ABC}                   tI10      I4/mmm              -         B2Ba1Rh2            #139      metal-boride; B2Ba1Rh2, ICSD #8155
A2B2C_tI10_139_d_e_a-001.{ABC}                   tI10      I4/mmm              -         B2Ir2Zn1            #139      metal-boride; B2Ir2Zn1, ICSD #71642
A2B2C_tI10_139_e_d_a-001.{ABC}                   tI10      I4/mmm              -         B2Co2Dy1            #139      metal-boride; B2Co2Dy1, ICSD #612907, equivalent to Cr2Si2Th (http://aflow.org/prototype-encyclopedia/A2B2C_tI10_139_d_e_a.html, part 3)
A2B3C7_hR12_166_c_bc_ah-001.{ABC}                hR12      R-3m                -         B2Ca3Ni7            #166      metal-boride; B2Ca3Ni7, ICSD #36505
A6B3C4_hR13_166_g_ac_2c-001.{ABC}                hR13      R-3m                -         B6Fe3Pr4            #166      metal-boride; B6Fe3Pr4, ICSD #614152
A4BC3_hP16_176_bh_c_h-001.{ABC}                  hP16      P6_{3}/m            -         B4Hf1Ir3            #176      metal-boride; B4Hf1Ir3, ICSD #1518
ABC2_hP12_180_d_c_i-001.{ABC}                    hP12      P6_{2}22            -         B1Ce1Pt2            #180      metal-boride; B1Ce1Pt2, ICSD #90403
ABC_hP9_189_ac_f_g-001.{ABC}                     hP9       P-62m               -         B1Fe1Nb1            #189      metal-boride; B1Fe1Nb1, ICSD #20298
A2B6C5_hP13_189_ab_i_df-001.{ABC}                hP13      P-62m               -         B2Rh6Sn5            #189      metal-boride; B2Rh6Sn5, ICSD #77352
ABC3_hP15_189_bc_f_fk-001.{ABC}                  hP15      P-62m               -         B1Li1Pt3            #189      metal-boride; B1Li1Pt3, ICSD #68091
ABC3_hP10_191_c_ab_i-001.{ABC}                   hP10      P6/mmm              -         B1Na1Pt3            #191      metal-boride; B1Na1Pt3, ICSD #68092
AB4C_hP12_191_d_ci_ab-001.{ABC}                  hP12      P6/mmm              -         B1Ni4Y1             #191      metal-boride; B1Ni4Y1, ICSD #16515
ABC4_hP12_191_c_ab_di-001.{ABC}                  hP12      P6/mmm              -         B1Ca1Ni4            #191      metal-boride; B1Ca1Ni4, ICSD #36504
A2B7C_cF40_225_c_bd_a.{ABC}                      cF40      Fm-3m               -         B2Pd7Y1             #225      metal-boride; B2Pd7Y1, ICSD #260147

*****************************************************************************************************************************************************
ANRL - metal oxides - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB_aP16_2_4i_4i-002.{AB}                         aP16      P-1                 -         Hg1O1               #2        metal-oxide; Hg1O1, ICSD #32561
AB2_mP12_4_2a_4a-001.{AB}                        mP12      P2_{1}              -         Mo1O2               #4        metal-oxide; Mo1O2, ICSD #36263
A4B_mC10_5_2c_a-001.{AB}                         mC10      C2                  -         O4Os1               #5        metal-oxide; O4Os1, ICSD #24672
A2B5_mC14_5_c_a2c-001.{AB}                       mC14      C2                  -         Nb2O5               #5        metal-oxide; Nb2O5, ICSD #51176
AB2_mC6_8_a_2a-001.{AB}                          mC6       Cm                  -         Co1O2               #8        metal-oxide; Co1O2, ICSD #95440
AB_mC8_9_a_a-001.{AB}                            mC8       Cc                  -         Cu1O1               #9        metal-oxide; Cu1O1, ICSD #69757
A2B_mP6_10_mn_ce-001.{AB}                        mP6       P2/m                -         O2V1                #10       metal-oxide; O2V1, ICSD #89471
A2B_mP12_10_2mo_im-001.{AB}                      mP12      P2/m                -         O2V1                #10       metal-oxide; O2V1, ICSD #1503
AB3_mP8_11_e_3e-001.{AB}                         mP8       P2_{1}/m            -         Mo1O3               #11       metal-oxide; Mo1O3, ICSD #80577
A3B2_mP10_11_3e_2e-001.{AB}                      mP10      P2_{1}/m            -         O3Pb2               #11       metal-oxide; O3Pb2, ICSD #36243
A2B_mP12_11_4e_2e-001.{AB}                       mP12      P2_{1}/m            -         O2Ti1               #11       metal-oxide; O2Ti1, ICSD #657748
A3B4_mP14_11_3e_4e-001.{AB}                      mP14      P2_{1}/m            -         O3Tl4               #11       metal-oxide; O3Tl4, ICSD #23478
A5B2_mP14_11_5e_2e-001.{AB}                      mP14      P2_{1}/m            -         O5V2                #11       metal-oxide; O5V2, ICSD #59960
AB2_mC6_12_a_i-002.{AB}                          mC6       C2/m                -         Hg1O2               #12       metal-oxide; Hg1O2, ICSD #48214
AB2_mC6_12_a_i-003.{AB}                          mC6       C2/m                -         Ni1O2               #12       metal-oxide; Ni1O2, ICSD #88720
AB2_mC12_12_i_2i-001.{AB}                        mC12      C2/m                -         Nd2O4               #12       metal-oxide; Nd2O4, ICSD #62228
AB2_mC12_12_i_2i-002.{AB}                        mC12      C2/m                -         O1Tl2               #12       metal-oxide; O1Tl2, ICSD #77699
AB2_mC18_12_ai_3i-001.{AB}                       mC18      C2/m                -         Mn1O2               #12       metal-oxide; Mn1O2, ICSD #150462
A3B7_mC20_12_di_a3i-001.{AB}                     mC20      C2/m                -         O3V7                #12       metal-oxide; O3V7, ICSD #77706
A2B3_mC20_12_2i_3i-002.{AB}                      mC20      C2/m                -         Nb2O3               #12       metal-oxide; Nb2O3, ICSD #263119
A3B7_mC20_12_ai_d3i-001.{AB}                     mC20      C2/m                -         O6V14               #12       metal-oxide; O6V14, ICSD #163494
A3B8_mC22_12_ai_2ij-001.{AB}                     mC22      C2/m                -         Cr3O8               #12       metal-oxide; Cr3O8, ICSD #155847
A2B_mC24_12_4i_2i-001.{AB}                       mC24      C2/m                -         O2V1                #12       metal-oxide; O2V1, ICSD #199
A2B_mC24_12_4i_2i-002.{AB}                       mC24      C2/m                -         O2Ti1               #12       metal-oxide; O2Ti1, ICSD #57154
A2B_mC24_12_4i_2i-003.{AB}                       mC24      C2/m                -         O2V1                #12       metal-oxide; O2V1, ICSD #73855
A2B_mC24_12_2ij_gi-001.{AB}                      mC24      C2/m                -         O2V1                #12       metal-oxide; O2V1, ICSD #34416
A5B8_mC26_12_ahi_2ij-001.{AB}                    mC26      C2/m                -         Mn5O8               #12       metal-oxide; Mn5O8, ICSD #16956
A5B3_mC32_12_5i_3i-002.{AB}                      mC32      C2/m                -         O5Ti3               #12       metal-oxide; O5Ti3, ICSD #26492
A4B3_mP14_13_2g_ag-001.{AB}                      mP14      P2/c                -         O4Sn3               #13       metal-oxide; O4Sn3, ICSD #174299
AB_mP8_14_ad_e-001.{AB}                          mP8       P2_{1}/c            -         Ag1O1               #14       metal-oxide; Ag1O1, ICSD #43741
A2B_mP12_14_2e_e-005.{AB}                        mP12      P2_{1}/c            -         O2Zr1               #14       metal-oxide; O2Zr1, ICSD #15983
A2B_mP12_14_2e_e-006.{AB}                        mP12      P2_{1}/c            -         O2Ti1               #14       metal-oxide; O2Ti1, ICSD #154035
A2B_mP12_14_2e_e-007.{AB}                        mP12      P2_{1}/c            -         O2Ti1               #14       metal-oxide; O2Ti1, ICSD #154036
A2B_mP12_14_2e_e-008.{AB}                        mP12      P2_{1}/c            -         O2Tc1               #14       metal-oxide; O2Tc1, ICSD #173153
A2B_mP12_14_2e_e-009.{AB}                        mP12      P2_{1}/c            -         O2Zr1               #14       metal-oxide; O2Zr1, ICSD #659226
A3B4_mP14_14_ae_2e-001.{AB}                      mP14      P2_{1}/c            -         Ag3O4               #14       metal-oxide; Ag3O4, ICSD #59225
A3B_mP16_14_3e_e-002.{AB}                        mP16      P2_{1}/c            -         O3Rb1               #14       metal-oxide; O3Rb1, ICSD #6094
A3B_mP16_14_3e_e-003.{AB}                        mP16      P2_{1}/c            -         O3Rb1               #14       metal-oxide; O3Rb1, ICSD #47164
A3B_mP16_14_3e_e-004.{AB}                        mP16      P2_{1}/c            -         O3Rb1               #14       metal-oxide; O3Rb1, ICSD #180568
A3B_mP16_14_3e_e-005.{AB}                        mP16      P2_{1}/c            -         O3W1                #14       metal-oxide; O3W1, ICSD #647640
AB_mC8_15_c_e-002.{AB}                           mC8       C2/c                -         Ag1O1               #15       metal-oxide; Ag1O1, ICSD #27659
AB_mC8_15_a_e-002.{AB}                           mC8       C2/c                -         Cu1O1               #15       metal-oxide; Cu1O1, ICSD #67850
AB2_mC12_15_c_f-001.{AB}                         mC12      C2/c                -         K1O2                #15       metal-oxide; K1O2, ICSD #38230
AB2_mC24_15_ae_2f-001.{AB}                       mC24      C2/c                -         Bi2O4               #15       metal-oxide; Bi2O4, ICSD #79500
A5B3_mC32_15_e2f_cf-001.{AB}                     mC32      C2/c                -         O5V3                #15       metal-oxide; O5V3, ICSD #15899
A5B3_mC32_15_e2f_af-001.{AB}                     mC32      C2/c                -         O5V3                #15       metal-oxide; O5V3, ICSD #32587
A8B3_oC22_21_achl_bg-001.{AB}                    oC22      C222                -         O8W3                #21       metal-oxide; O8W3, ICSD #73719
AB_oP8_29_a_a-001.{AB}                           oP8       Pca2_{1}            -         O1Pb1               #29       metal-oxide; O1Pb1, ICSD #36250
A2B_oP12_29_2a_a-002.{AB}                        oP12      Pca2_{1}            -         O2Zr1               #29       metal-oxide; O2Zr1, ICSD #67004
AB_oP4_31_a_a-001.{AB}                           oP4       Pmn2_{1}            -         O1Sn1               #31       metal-oxide; O1Sn1, ICSD #20624
A5B2_oP14_31_a2b_b-001.{AB}                      oP14      Pmn2_{1}            D8_{7}    O5V2                #31       metal-oxide; O5V2, ICSD #15984, equivalent to O5V2 (shcherbnaite (obsolete), http://aflow.org/prototype-encyclopedia/A5B2_oP14_31_a2b_b.html, part 3)
A5B2_oP14_31_a2b_b-002.{AB}                      oP14      Pmn2_{1}            -         O5V2                #31       metal-oxide; O5V2, ICSD #29140
A5B2_oP14_31_a2b_b-003.{AB}                      oP14      Pmn2_{1}            -         O5V2                #31       metal-oxide; O5V2, ICSD #41030
A2B_oC24_35_abef_de-001.{AB}                     oC24      Cmm2                -         O2Ti1               #35       metal-oxide; O2Ti1, ICSD #97008
AB_oC16_36_2a_2a-002.{AB}                        oC16      Cmc2_{1}            -         O1Sn1               #36       metal-oxide; O1Sn1, ICSD #424729
AB_oC16_36_b_2a-001.{AB}                         oC16      Cmc2_{1}            -         O1Sn1               #36       metal-oxide; O1Sn1, ICSD #60619
A8B_oC18_38_2cde_a-001.{AB}                      oC18      Amm2                -         Cu8O1               #38       metal-oxide; Cu8O1, ICSD #62764
AB3_oC32_40_c_3c-001.{AB}                        oC32      Ama2                -         Cr1O3               #40       metal-oxide; Cr1O3, ICSD #16031
A2B3_oF40_43_b_ab-001.{AB}                       oF40      Fdd2                -         Au2O3               #43       metal-oxide; Au2O3, ICSD #8014, equivalent to Ag2O3 (http://aflow.org/prototype-encyclopedia/A2B3_oF40_43_b_ab.html, part 3)
AB3_oI8_44_a_ac-001.{AB}                         oI8       Imm2                -         Li1O3               #44       metal-oxide; Li1O3, ICSD #180565
AB_oP2_47_a_e-001.AB                             oP2       Pmmm                -         Nb1O1               #47       metal-oxide; Nb1O1, ICSD #95729
AB2_oP3_47_a_m-001.{AB}                          oP3       Pmmm                -         Cu1O2               #47       metal-oxide; Cu1O2, ICSD #54126
AB4_oP5_47_e_bdr-001.{AB}                        oP5       Pmmm                -         O1Ta4               #47       metal-oxide; O1Ta4, ICSD #76022
A5B2_oP7_47_behq_ad-001.{AB}                     oP7       Pmmm                -         O5Ta2               #47       metal-oxide; O5Ta2, ICSD #95462
A4B3_oP14_55_gh_bg-001.{AB}                      oP14      Pbam                -         O4Pb3               #55       metal-oxide; O4Pb3, ICSD #97282
AB2_oP12_59_e_2e-001.{AB}                        oP12      Pmmn                -         Mn1O2               #59       metal-oxide; Mn1O2, ICSD #54114
AB2_oP12_59_e_2e-002.{AB}                        oP12      Pmmn                -         Hf1O2               #59       metal-oxide; Hf1O2, ICSD #83863
A5B2_oP14_59_a2e_e-002.{AB}                      oP14      Pmmn                -         O5V2                #59       metal-oxide; O5V2, ICSD #43132
A2B_oP12_60_d_c-002.{AB}                         oP12      Pbcn                -         O2Ti1               #60       metal-oxide; O2Ti1, ICSD #15328
A2B_oP12_60_d_c-003.{AB}                         oP12      Pbcn                -         O2Pb1               #60       metal-oxide; O2Pb1, ICSD #20362, equivalent to alpha-PbO2 (http://aflow.org/prototype-encyclopedia/A2B_oP12_60_d_c.html, part 3)
A2B_oP12_60_d_c-004.{AB}                         oP12      Pbcn                -         O2Re1               #60       metal-oxide; O2Re1, ICSD #24060
AB2_oP12_61_a_c-001.{AB}                         oP12      Pbca                -         Hg1O2               #61       metal-oxide; Hg1O2, ICSD #24774
A2B_oP12_62_2c_c-006.{AB}                        oP12      Pnma                -         O2Zr1               #62       metal-oxide; O2Zr1, ICSD #56696 equivalent to 652.AB,652.BA [HTQC]
AB2_oP12_62_c_2c-001.{AB}                        oP12      Pnma                -         Mn1O2               #62       metal-oxide; Mn1O2, ICSD #171866
A2B_oP12_62_2c_c-007.{AB}                        oP12      Pnma                -         O2Sn1               #62       metal-oxide; O2Sn1, ICSD #181282
A2B_oP12_62_2c_c-008.{AB}                        oP12      Pnma                -         O2Ti1               #62       metal-oxide; O2Ti1, ICSD #182577
A2B_oP12_62_2c_c-009.{AB}                        oP12      Pnma                -         O2Pr1               #62       metal-oxide; O2Pr1, ICSD #380398
AB3_oP16_62_c_3c-002.{AB}                        oP16      Pnma                -         Mo1O3               #62       metal-oxide; Mo1O3, ICSD #644065 equivalent to 698.AB,698.BA [HTQC]
A2B3_oC20_63_ac_cf-001.{AB}                      oC20      Cmcm                -         Al2O3               #63       metal-oxide; Al2O3, ICSD #161062
A2B3_oC20_63_ac_cf-002.{AB}                      oC20      Cmcm                -         Ga2O3               #63       metal-oxide; Ga2O3, ICSD #162252 equivalent to 1274.AB,1274.BA [HTQC]
A3B4_oC28_63_cf_acf-001.{AB}                     oC28      Cmcm                -         Fe3O4               #63       metal-oxide; Fe3O4, ICSD #263010
AB2_oC6_65_a_j-001.{AB}                          oC6       Cmmm                -         Ba1O2               #65       metal-oxide; Ba1O2, ICSD #180398
AB2_oF12_69_a_h-001.{AB}                         oF12      Fmmm                -         Cu1O2               #69       metal-oxide; Cu1O2, ICSD #85080
AB2_oF12_69_a_g-001.{AB}                         oF12      Fmmm                -         Cu1O2               #69       metal-oxide; Cu1O2, ICSD #150886
AB_oI8_71_g_e-001.{AB}                           oI8       Immm                -         O2Rb2               #71       metal-oxide; O2Rb2, ICSD #25528, equivalent to CsO (http://aflow.org/prototype-encyclopedia/AB_oI8_71_g_i.html, part 3)
AB_oI32_72_2j_2j-001.AB                          oI32      Ibam                -         Mg1O1               #72       metal-oxide; Mg1O1, ICSD #181462
A3B4_oI28_74_ace_hi-001.{AB}                     oI28      Imma                -         Fe3O4               #74       metal-oxide; Fe3O4, ICSD #31156
AB2_tI24_87_h_2h-001.{AB}                        tI24      I4/m                -         Mn1O2               #87       metal-oxide; Mn1O2, ICSD #20227
AB_tI32_88_cd_f-001.{AB}                         tI32      I4_{1}/a            -         Ag1O1               #88       metal-oxide; Ag1O1, ICSD #202055
A3B_tP16_113_ef_e-001.{AB}                       tP16      P-42_{1}m           -         O3W1                #113      metal-oxide; O3W1, ICSD #86144
A2B3_tP5_115_g_abc-001.{AB}                      tP5       P-4m2               -         Bi2O3               #115      metal-oxide; Bi2O3, ICSD #168808
A2B_tP3_123_g_a-002.{AB}                         tP3       P4/mmm              -         O2Zr1               #123      metal-oxide; O2Zr1, ICSD #92091
A3B_tP8_127_bg_a-001.{AB}                        tP8       P4/mbm              -         O3Re1               #127      metal-oxide; O3Re1, ICSD #77680
AB_tP4_129_c_a-001.{AB}                          tP4       P4/nmm              -         Ba1O1               #129      metal-oxide; Ba1O1, ICSD #15301
AB2_tP6_129_c_ac-001.{AB}                        tP6       P4/nmm              -         Gd1O2               #129      metal-oxide; Gd1O2, ICSD #6318
AB2_tP6_129_c_ac-002.{AB}                        tP6       P4/nmm              -         Eu1O2               #129      metal-oxide; Eu1O2, ICSD #631464
A3B_tP8_129_cd_c-001.{AB}                        tP8       P4/nmm              -         O3W1                #129      metal-oxide; O3W1, ICSD #27961
A3B_tP16_130_cf_c-001.{AB}                       tP16      P4/ncc              -         O3W1                #130      metal-oxide; O3W1, ICSD #50732, equivalent to alpha-WO3 (http://aflow.org/prototype-encyclopedia/A3B_tP16_130_cf_c.html, part 3)
A2B3_tP10_132_i_be-001.{AB}                      tP10      P4_{2}/mcm          -         Bi2O3               #132      metal-oxide; Bi2O3, ICSD #168807
A2B_tP6_136_f_a-002.{AB}                         tP6       P4_{2}/mnm          -         O2W1                #136      metal-oxide; O2W1, ICSD #647647
A2B_tP6_137_d_a-002.{AB}                         tP6       P4_{2}/nmc          -         O2Zr1               #137      metal-oxide; O2Zr1, ICSD #93028
AB_tI4_139_b_a-001.AB                            tI4       I4/mmm              -         O1Pd1               #139      metal-oxide; O1Pd1, ICSD #41617
AB_tI4_139_a_b-002.AB                            tI4       I4/mmm              -/L2_{0}  Co1O1/FeC_{x}       #139      metal-oxide; Co1O1, ICSD #174027, equivalent to "Martensite Type" FeC_{x} (L2_{0}, http://aflow.org/prototype-encyclopedia/AB_tI4_139_b_a.html, part 3)
AB2_tI6_139_a_e-005.{AB}                         tI6       I4/mmm              -         Ca1O2               #139      metal-oxide; Ca1O2, ICSD #20275
AB2_tI6_139_a_e-006.{AB}                         tI6       I4/mmm              -         Ba1O2               #139      metal-oxide; Ba1O2, ICSD #24729
AB3_tI32_140_ab_hk-001.{AB}                      tI32      I4/mcm              -         K1O3                #140      metal-oxide; K1O3, ICSD #28347
A2B_tI12_141_e_a-004.{AB}                        tI12      I4_{1}/amd          -         O2Ti1               #141      metal-oxide; O2Ti1, ICSD #200392
A2B_tI24_141_h_c-001.{AB}                        tI24      I4_{1}/amd          -         O2Ti1               #141      metal-oxide; O2Ti1, ICSD #93098
A4B3_tI28_141_cd_ae-001.{AB}                     tI28      I4_{1}/amd          -         Cu4O3               #141      metal-oxide; Cu4O3, ICSD #77675
A3B4_tI28_141_ad_h-002.{AB}                      tI28      I4_{1}/amd          -         Mn3O4               #141      metal-oxide; Mn3O4, ICSD #643198
A6B_hR7_148_f_a-001.{AB}                         hR7       R-3                 -         Hf1O0.1667          #148      metal-oxide; Hf1O0.1667, ICSD #174039
A2B3_hP5_150_d_e-001.{AB}                        hP5       P321                -         La2O3               #150      metal-oxide; La2O3, ICSD #26864
AB_hP6_152_b_a-001.{AB}                          hP6       P3_{1}21            -         Hg1O1               #152      metal-oxide; Hg1O1, ICSD #24062
A2B_hP9_152_c_a-002.{AB}                         hP9       P3_{1}21            -         O2Ti1               #152      metal-oxide; O2Ti1, ICSD #41493
AB3_hR16_155_2c_def-001.{AB}                     hR16      R32                 -         O1Zr3               #155      metal-oxide; O1Zr3, ICSD #23402
AB_hR2_160_a_a-001.AB                            hR2       R3m                 -         Bi1O1               #160      metal-oxide; Bi1O1, ICSD #30361
A2B3_hR10_160_ab_2b-001.{AB}                     hR10      R3m                 -         Bi2O3               #160      metal-oxide; Bi2O3, ICSD #168810
AB6_hP7_162_a_k-001.{AB}                         hP7       P-31m               -         O1Ti6               #162      metal-oxide; O1Ti6, ICSD #20042
A3B_hP8_162_k_d-001.{AB}                         hP8       P-31m               -         Ag6O2               #162      metal-oxide; Ag6O2, ICSD #26557
AB3_hP8_162_ab_k-001.{AB}                        hP8       P-31m               -         O1Ti3               #162      metal-oxide; O1Ti3, ICSD #23575
AB6_hP14_163_c_i-001.{AB}                        hP14      P-31c               -         O1Ti6               #163      metal-oxide; O1Ti6, ICSD #17009
AB3_hP16_163_ac_i-001.{AB}                       hP16      P-31c               -         O1Ti3               #163      metal-oxide; O1Ti3, ICSD #20041
A2B_hP3_164_d_a-004.{AB}                         hP3       P-3m1               -         O2Pt1               #164      metal-oxide; O2Pt1, ICSD #76431
A3B2_hP5_164_ad_d-005.{AB}                       hP5       P-3m1               -         O3Sc2               #164      metal-oxide; O3Sc2, ICSD #160203
AB2_hP6_164_d_abd-001.{AB}                       hP6       P-3m1               -         Ce1O2               #164      metal-oxide; Ce1O2, ICSD #189287
AB_hR2_166_a_b-002.AB                            hR2       R-3m                -         Fe1O1               #166      metal-oxide; Fe1O1, ICSD #82236
A2B_hR3_166_c_a-004.{AB}                         hR3       R-3m                -         Li2O1               #166      metal-oxide; Li2O1, ICSD #108886
A3B4_hR14_166_acd_ch-001.{AB}                    hR14      R-3m                -         Fe3O4               #166      metal-oxide; Fe3O4, ICSD #92356
AB_hR16_166_ab3c_4c-001.{AB}                     hR16      R-3m                -         Mg1O1               #166      metal-oxide; Mg1O1, ICSD #181459
AB3_hR8_167_b_e-001.{AB}                         hR8       R-3c                -         O1Zr3               #167      metal-oxide; O1Zr3, ICSD #27023 equivalent to 1265.AB,1265.BA [HTQC]
A2B3_hR10_167_c_e-003.{AB}                       hR10      R-3c                -         Al2O3               #167      metal-oxide; Al2O3, ICSD #43732
A2B3_hR10_167_c_e-004.{AB}                       hR10      R-3c                -         Al2O3               #167      metal-oxide; Al2O3, ICSD #99783
A2B_hP6_186_ab_b-001.{AB}                        hP6       P6_{3}mc            -         O2Pt1               #186      metal-oxide; O2Pt1, ICSD #24923
AB_hP2_187_d_a-003.AB                            hP2       P-6m2               -         Mg1O1               #187      metal-oxide; Mg1O1, ICSD #166273
A3B_hP12_191_gl_f-001.{AB}                       hP12      P6/mmm              -         O3W1                #191      metal-oxide; O3W1, ICSD #32001 equivalent to 699.AB,699.BA [HTQC], equivalent to hexagonal O3W (http://aflow.org/prototype-encyclopedia/A3B_hP12_191_gl_f.html, part 3)
A3B_hP8_193_g_b-002.{AB}                         hP8       P6_{3}/mcm          -         Cs3O1               #193      metal-oxide; Cs3O1, ICSD #15695
AB_hP8_194_f_f-001.{AB}                          hP8       P6_{3}/mmc          -         Li2O2               #194      metal-oxide; Li2O2, ICSD #24143
AB_hP8_194_f_f-002.{AB}                          hP8       P6_{3}/mmc          -         Li2O2               #194      metal-oxide; Li2O2, ICSD #50658
A3B_cI32_204_g_c-002.{AB}                        cI32      Im-3                -         O3Re1               #204      metal-oxide; O3Re1, ICSD #55462
A3B_cI32_204_g_c-003.{AB}                        cI32      Im-3                -         O3Re1               #204      metal-oxide; O3Re1, ICSD #55469
A3B4_cF56_216_abe_2e-001.{AB}                    cF56      F-43m               -         Fe3O4               #216      metal-oxide; Fe3O4, ICSD #65338
AB2_cI24_217_c_abc-001.{AB}                      cI24      I-43m               -         O1Ta2               #217      metal-oxide; O1Ta2, ICSD #28387
AB_cP12_223_c_d.AB                               cP12      Pm-3n               -         Mg1O1               #223      metal-oxide; Mg1O1, ICSD #181465
A4B3_cP14_223_e_c.{AB}                           cP14      Pm-3n               -         O4Pt3               #223      metal-oxide; O4Pt3, ICSD #30444
A2B3_cP10_224_b_d.{AB}                           cP10      Pn-3m               D5_{5}    Ag2O3 (Mg3P2)       #224      metal-oxide; Ag2O3, ICSD #15999, equivalent to Mg3P2 (D5_{5}, http://aflow.org/prototype-encyclopedia/A3B2_cP10_224_d_b.html, part 3)
A3B4_cF56_227_ad_e-002.{AB}                      cF56      Fd-3m               -         Co3O4               #227      metal-oxide; Co3O4, ICSD #63164

*****************************************************************************************************************************************************
ANRL - metal oxides - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB3C_aP5_1_a_3a_a-001.{ABC}                      aP5       P1                  -         Cu1O3V1             #1        metal-oxide; Cu1O3V1, ICSD #9414
ABC3_aP10_1_2a_2a_6a-001.{ABC}                   aP10      P1                  -         Ca1Fe1O3            #1        metal-oxide; Ca1Fe1O3, ICSD #92339
A2B4C_aP14_1_4a_8a_2a-001.{ABC}                  aP14      P1                  -         Al2O4Pb1            #1        metal-oxide; Al2O4Pb1, ICSD #33532
A2BC4_aP14_1_4a_2a_8a-001.{ABC}                  aP14      P1                  -         Cr2Ni1O4            #1        metal-oxide; Cr2Ni1O4, ICSD #280061
AB4C3_aP16_1_2a_8a_6a-001.{ABC}                  aP16      P1                  -         Co1Na4O3            #1        metal-oxide; Co1Na4O3, ICSD #10473
AB6C2_aP9_2_a_3i_i-001.{ABC}                     aP9       P-1                 -         Cu1O6V2             #2        metal-oxide; Cu1O6V2, ICSD #28151
AB3C_aP10_2_i_3i_i-001.ABC/BCA/CAB               aP10      P-1                 -         Mn1O3Sn1            #2        metal-oxide; Mn1O3Sn1, ICSD #29203
AB3C_aP10_2_i_3i_i-002.{ABC}                     aP10      P-1                 -         Hg1O3V1             #2        metal-oxide; Hg1O3V1, ICSD #82242
AB3C6_aP10_2_a_bi_3i-001.{ABC}                   aP10      P-1                 -         Au1Gd3O6            #2        metal-oxide; Au1Gd3O6, ICSD #411500
AB4C_aP12_2_i_4i_i-001.{ABC}                     aP12      P-1                 -         Cu1O4W1             #2        metal-oxide; Cu1O4W1, ICSD #4189
A3B8C2_aP13_2_ai_4i_i-001.{ABC}                  aP13      P-1                 -         Cu3O8V2             #2        metal-oxide; Cu3O8V2, ICSD #27184
AB3C3_aP14_2_i_3i_3i-001.{ABC}                   aP14      P-1                 -         Co2Na6O6            #2        metal-oxide; Co2Na6O6, ICSD #99580
A4BC2_aP14_2_4i_i_abcf-001.{ABC}                 aP14      P-1                 -         O4Pb1Pt2            #2        metal-oxide; O4Pb1Pt2, ICSD #59657
A4BC3_aP16_2_4i_i_3i-001.{ABC}                   aP16      P-1                 -         Cs4Fe1O3            #2        metal-oxide; Cs4Fe1O3, ICSD #423336
ABC2_mP16_4_2a_2a_4a-001.{ABC}                   mP16      P2_{1}              -         Ca1Cu1O2            #4        metal-oxide; Ca1Cu1O2, ICSD #84868
AB6C2_mC18_5_a_3c_c-001.{ABC}                    mC18      C2                  -         Cu1O6V2             #5        metal-oxide; Cu1O6V2, ICSD #21067
A6B2C_mC18_5_3c_c_a-001.{ABC}                    mC18      C2                  -         O6V2Zn1             #5        metal-oxide; O6V2Zn1, ICSD #26998
ABC4_mC24_5_c_c_4c-001.{ABC}                     mC24      C2                  -         Fe1Nb1O4            #5        metal-oxide; Fe1Nb1O4, ICSD #429
A8B3C2_mC26_5_4c_ac_c-001.{ABC}                  mC26      C2                  -         O8Pb3V2             #5        metal-oxide; O8Pb3V2, ICSD #29359
AB4C2_mC28_5_c_4c_abc-001.{ABC}                  mC28      C2                  -         Mo1O4Tl2            #5        metal-oxide; Mo1O4Tl2, ICSD #280608
AB3C_mP5_6_b_2ab_a-001.{ABC}                     mP5       Pm                  -         Ba1O3Ti1            #6        metal-oxide; Ba1O3Ti1, ICSD #186460
A2BC3_mC12_8_2a_a_3a-001.{ABC}                   mC12      Cm                  -         Ca2Co1O3            #8        metal-oxide; Ca2Co1O3, ICSD #95439
AB5C5_mC22_8_a_5a_5a-001.{ABC}                   mC22      Cm                  -         Bi1Li5O5            #8        metal-oxide; Bi1Li5O5, ICSD #203031
ABC2_mC24_8_ab_ab_2a2b-001.ABC/BCA/CAB           mC24      Cm                  -         Li1Ni1O2            #8        metal-oxide; Li1Ni1O2, ICSD #164853
A2B6C5_mC26_8_2a_2a2b_a2b-001.{ABC}              mC26      Cm                  -         Fe2K6O5             #8        metal-oxide; Fe2K6O5, ICSD #174312
A6B2C5_mC26_8_2a2b_2a_a2b-001.{ABC}              mC26      Cm                  -         Cs6Fe2O5            #8        metal-oxide; Cs6Fe2O5, ICSD #174313
ABC3_mC20_9_a_a_3a-001.{ABC}                     mC20      Cc                  -         Bi1Fe1O3            #9        metal-oxide; Bi1Fe1O3, ICSD #247765
A2B3C_mC24_9_2a_3a_a-001.{ABC}                   mC24      Cc                  -         Li2O3Zr1            #9        metal-oxide; Li2O3Zr1, ICSD #31941
AB4C3_mC32_9_a_4a_3a-001.{ABC}                   mC32      Cc                  -         Fe1Na4O3            #9        metal-oxide; Fe1Na4O3, ICSD #1410
AB4C3_mC32_9_a_4a_3a-002.{ABC}                   mC32      Cc                  -         Co1Na4O3            #9        metal-oxide; Co1Na4O3, ICSD #10501
A4B3C_mC32_9_4a_3a_a-001.{ABC}                   mC32      Cc                  -         Na4O3Sn1            #9        metal-oxide; Na4O3Sn1, ICSD #49624
ABC3_mP10_11_a_e_d2e-001.{ABC}                   mP10      P2_{1}/m            -         Bi1Fe1O3            #11       metal-oxide; Bi1Fe1O3, ICSD #162834
AB2C2_mP10_11_e_2e_2e-002.{ABC}                  mP10      P2_{1}/m            -         Cu1Na2O2            #11       metal-oxide; Cu1Na2O2, ICSD #422751
AB4C_mP12_11_e_aef_e-001.{ABC}                   mP12      P2_{1}/m            -         Ba1O4W1             #11       metal-oxide; Ba1O4W1, ICSD #155516
A2BC4_mP14_11_2e_e_4e-001.{ABC}                  mP14      P2_{1}/m            -         Al2Ca1O4            #11       metal-oxide; Al2Ca1O4, ICSD #172780
ABC2_mC8_12_a_c_i-001.{ABC}                      mC8       C2/m                -         Cu1Na1O2            #12       metal-oxide; Cu1Na1O2, ICSD #15098
ABC2_mC8_12_a_c_i-002.{ABC}                      mC8       C2/m                -         Mn1Na1O2            #12       metal-oxide; Mn1Na1O2, ICSD #16270
ABC2_mC8_12_a_c_i-003.{ABC}                      mC8       C2/m                -         Cu1Li1O2            #12       metal-oxide; Cu1Li1O2, ICSD #74978
ABC2_mC8_12_c_a_i-001.{ABC}                      mC8       C2/m                -         Na1Ni1O2            #12       metal-oxide; Na1Ni1O2, ICSD #26609
ABC2_mC8_12_c_a_i-002.{ABC}                      mC8       C2/m                -         Cu1Mn1O2            #12       metal-oxide; Cu1Mn1O2, ICSD #30379
ABC2_mC8_12_a_d_i-001.{ABC}                      mC8       C2/m                -         Cd1Hg1O2            #12       metal-oxide; Cd1Hg1O2, ICSD #74848
ABC2_mC8_12_a_d_i-002.{ABC}                      mC8       C2/m                -         Ag1Cu1O2            #12       metal-oxide; Ag1Cu1O2, ICSD #95089
AB2C_mC8_12_a_i_b-001.{ABC}                      mC8       C2/m                -         Na1O2V1             #12       metal-oxide; Na1O2V1, ICSD #420136
AB2C_mC8_12_a_i_c-001.{ABC}                      mC8       C2/m                -         Na1O2V1             #12       metal-oxide; Na1O2V1, ICSD #420138
A2BC2_mC10_12_i_a_i-002.{ABC}                    mC10      C2/m                -         Ag2Ni1O2            #12       metal-oxide; Ag2Ni1O2, ICSD #172559
AB2C2_mC10_12_a_i_i-001.{ABC}                    mC10      C2/m                -         Cu1Li2O2            #12       metal-oxide; Cu1Li2O2, ICSD #174134
ABC_mC12_12_i_i_i-002.ABC/BCA/CAB                mC12      C2/m                -         K1O1Tl1             #12       metal-oxide; K1O1Tl1, ICSD #1570
ABC2_mC16_12_ad_i_2i-001.{ABC}                   mC16      C2/m                -         Li1Mo1O2            #12       metal-oxide; Li1Mo1O2, ICSD #165326
AB4C4_mC18_12_a_2i_j-001.{ABC}                   mC18      C2/m                -         Ir1K4O4             #12       metal-oxide; Ir1K4O4, ICSD #47223
A2B3C4_mC18_12_i_ai_2i-001.{ABC}                 mC18      C2/m                -         Cu2Li3O4            #12       metal-oxide; Cu2Li3O4, ICSD #66509
A2B5C2_mC18_12_i_a2i_i-001.{ABC}                 mC18      C2/m                -         K2O5Ti2             #12       metal-oxide; K2O5Ti2, ICSD #36097, equivalent to K2O5Ti2 (http://aflow.org/prototype-encyclopedia/A2B5C2_mC18_12_i_a2i_i.html, part 3)
A4BC4_mC18_12_2i_a_j-001.{ABC}                   mC18      C2/m                -         Cs4Ir1O4            #12       metal-oxide; Cs4Ir1O4, ICSD #72287
AB6C2_mC18_12_a_3i_i-003.{ABC}                   mC18      C2/m                -         Co1O6V2             #12       metal-oxide; Co1O6V2, ICSD #188911
AB3C_mC20_12_i_ghi_e-001.{ABC}                   mC20      C2/m                -         Ba1O3Pb1            #12       metal-oxide; Ba1O3Pb1, ICSD #94313
AB8C2_mC22_12_b_ac3i_i-001.{ABC}                 mC22      C2/m                -         Mo1O8V2             #12       metal-oxide; Mo1O8V2, ICSD #28471
A2B8C_mC22_12_i_2ij_a-001.{ABC}                  mC22      C2/m                -         Mo2O8Zr1            #12       metal-oxide; Mo2O8Zr1, ICSD #280433
A7B2C2_mC22_12_aij_i_h-001.{ABC}                 mC22      C2/m                -         O7V2Zn2             #12       metal-oxide; O7V2Zn2, ICSD #250002
A3BC8_mC24_12_ci_a_2ij-001.{ABC}                 mC24      C2/m                -         Cr3K1O8             #12       metal-oxide; Cr3K1O8, ICSD #15898
A3BC8_mC24_12_ci_a_2ij-002.{ABC}                 mC24      C2/m                -         Cr3Na1O8            #12       metal-oxide; Cr3Na1O8, ICSD #82620
A3BC8_mC24_12_ci_a_2ij-003.{ABC}                 mC24      C2/m                -         Cr3In1O8            #12       metal-oxide; Cr3In1O8, ICSD #155504
A3BC8_mC24_12_ci_a_2ij-004.{ABC}                 mC24      C2/m                -         Cr3K1O8             #12       metal-oxide; Cr3K1O8, ICSD #246524
ABC4_mC24_12_i_i_4i-001.{ABC}                    mC24      C2/m                -         Fe1Nb1O4            #12       metal-oxide; Fe1Nb1O4, ICSD #14016, equivalent to AlNbO4 (http://aflow.org/prototype-encyclopedia/ABC4_mC24_12_i_i_4i.html, part 3)
AB4C_mC24_12_i_2ij_g-001.{ABC}                   mC24      C2/m                -         Al1O4W1             #12       metal-oxide; Al1O4W1, ICSD #4164
AB4C_mC24_12_i_2ij_g-002.{ABC}                   mC24      C2/m                -         Cr1O4Re1            #12       metal-oxide; Cr1O4Re1, ICSD #173875
A2BC3_mC24_12_acg_h_ij-001.{ABC}                 mC24      C2/m                -         Li2Ni1O3            #12       metal-oxide; Li2Ni1O3, ICSD #153094
A4B7C_mC24_12_2i_b3i_a-001.{ABC}                 mC24      C2/m                -         La4O7Pd1            #12       metal-oxide; La4O7Pd1, ICSD #65032
A5B6C_mC24_12_dgh_ij_a-001.{ABC}                 mC24      C2/m                -         Li5O6Re1            #12       metal-oxide; Li5O6Re1, ICSD #38381
A5B6C_mC24_12_dgh_ij_a-002.{ABC}                 mC24      C2/m                -         Na5O6Re1            #12       metal-oxide; Na5O6Re1, ICSD #38382
A3B8C_mC24_12_ci_2ij_a-001.{ABC}                 mC24      C2/m                -         Cr3O8Tl1            #12       metal-oxide; Cr3O8Tl1, ICSD #155505
AB3C8_mC24_12_a_ci_2ij-001.{ABC}                 mC24      C2/m                -         Ag1Cr3O8            #12       metal-oxide; Ag1Cr3O8, ICSD #155507
A3BC8_mC24_12_di_a_2ij-001.{ABC}                 mC24      C2/m                -         Cr3Li1O8            #12       metal-oxide; Cr3Li1O8, ICSD #155851
A4BC8_mC26_12_2i_a_4i-001.{ABC}                  mC26      C2/m                -         Cr8K2O16            #12       metal-oxide; Cr8K2O16, ICSD #491
A2B3C8_mC26_12_i_ah_2ij-001.{ABC}                mC26      C2/m                -         Cu2Mn3O8            #12       metal-oxide; Cu2Mn3O8, ICSD #971
A2B3C8_mC26_12_i_ah_2ij-002.{ABC}                mC26      C2/m                -         Cd2Mn3O8            #12       metal-oxide; Cd2Mn3O8, ICSD #16957
A3B8C2_mC26_12_ai_2ij_i-001.{ABC}                mC26      C2/m                -         Ca3O8V2             #12       metal-oxide; Ca3O8V2, ICSD #412273
A6B2C5_mC26_12_ghi_i_aj-001.{ABC}                mC26      C2/m                -         Cs6Fe2O5            #12       metal-oxide; Cs6Fe2O5, ICSD #73134
A4B3C6_mC26_12_2i_ai_hj-001.{ABC}                mC26      C2/m                -         K4Ni3O6             #12       metal-oxide; K4Ni3O6, ICSD #426450
A4B7C2_mC26_12_aci_bij_i-001.{ABC}               mC26      C2/m                -         K4O7V2              #12       metal-oxide; K4O7V2, ICSD #250388
AB4C8_mC26_12_a_2i_4i-001.{ABC}                  mC26      C2/m                -         Ba2Mn8O16           #12       metal-oxide; Ba2Mn8O16, ICSD #62096
A4BC8_mC26_12_2i_a_4i-002.{ABC}                  mC26      C2/m                -         Ir1K0.25O2          #12       metal-oxide; Ir1K0.25O2, ICSD #80336
A4B8C_mC26_12_2i_4i_a-001.{ABC}                  mC26      C2/m                -         Ga4O8Ti1            #12       metal-oxide; Ga4O8Ti1, ICSD #155638
AB3C3_mC28_12_i_ghi_ij-001.{ABC}                 mC28      C2/m                -         Ga2K6O6             #12       metal-oxide; Ga2K6O6, ICSD #2269
AB2C4_mC28_12_i_2i_4i-001.{ABC}                  mC28      C2/m                -         Hg1Ho2O4            #12       metal-oxide; Hg1Ho2O4, ICSD #69734
A2BC4_mC28_12_2i_i_2ij-001.{ABC}                 mC28      C2/m                -         K2Mo1O4             #12       metal-oxide; K2Mo1O4, ICSD #16154
A2B4C_mC28_12_2i_4i_i-001.{ABC}                  mC28      C2/m                -         Bi2O4Sr1            #12       metal-oxide; Bi2O4Sr1, ICSD #80668
A3B10C2_mC30_12_ah_3ij_i-001.{ABC}               mC30      C2/m                -         La3O10Os2           #12       metal-oxide; La3O10Os2, ICSD #10104
A10B3C2_mC30_12_3ij_ah_i-001.{ABC}               mC30      C2/m                -         O10Ru3Sr2           #12       metal-oxide; O10Ru3Sr2, ICSD #50707
AB5C2_mC32_12_i_g2ij_2i-001.{ABC}                mC32      C2/m                -         Cr1O5Pb2            #12       metal-oxide; Cr1O5Pb2, ICSD #29269
AB5C2_mC32_12_i_g2ij_2i-002.{ABC}                mC32      C2/m                -         Cr1O5Pb2            #12       metal-oxide; Cr1O5Pb2, ICSD #201562
A5B2C_mC32_12_g2ij_2i_i-001.{ABC}                mC32      C2/m                -         O5Pb2W1             #12       metal-oxide; O5Pb2W1, ICSD #61399
AB4C_mP12_13_f_2g_e-003.{ABC}                    mP12      P2/c                -         Mg1O4W1             #13       metal-oxide; Mg1O4W1, ICSD #36310
AB4C_mP12_13_f_2g_e-004.{ABC}                    mP12      P2/c                -         Cu1O4W1             #13       metal-oxide; Cu1O4W1, ICSD #182751
AB4C_mP12_13_e_2g_f-001.{ABC}                    mP12      P2/c                -         Mo1O4Zn1            #13       metal-oxide; Mo1O4Zn1, ICSD #236416
AB4C_mP12_13_e_2g_f-002.{ABC}                    mP12      P2/c                -         Nd1O4Ta1            #13       metal-oxide; Nd1O4Ta1, ICSD #415430
AB3C4_mP16_13_e_e2f_2g-001.{ABC}                 mP16      P2/c                -         Bi1Na3O4            #13       metal-oxide; Bi1Na3O4, ICSD #10319
A2B2C3_mP14_14_e_e_ae-001.{ABC}                  mP14      P2_{1}/c            -         Cd2K2O3             #14       metal-oxide; Cd2K2O3, ICSD #16223
AB2C_mP16_14_e_2e_e-001.{ABC}                    mP16      P2_{1}/c            -         Li1O2Y1             #14       metal-oxide; Li1O2Y1, ICSD #45511
ABC2_mP16_14_e_e_2e-001.{ABC}                    mP16      P2_{1}/c            -         Ce1Li1O2            #14       metal-oxide; Ce1Li1O2, ICSD #47116
ABC2_mC16_15_e_e_f-001.{ABC}                     mC16      C2/c                -         Bi1Na1O2            #15       metal-oxide; Bi1Na1O2, ICSD #10317
ABC2_mC16_15_e_e_f-002.{ABC}                     mC16      C2/c                -         Bi1K1O2             #15       metal-oxide; Bi1K1O2, ICSD #16239
ABC2_mC16_15_e_e_f-003.{ABC}                     mC16      C2/c                -         Bi1Cs1O2            #15       metal-oxide; Bi1Cs1O2, ICSD #406564
AB2C_mC16_15_e_f_e-001.{ABC}                     mC16      C2/c                -         Bi1O2Rb1            #15       metal-oxide; Bi1O2Rb1, ICSD #407208
ABC3_mC20_15_a_e_ef-001.{ABC}                    mC20      C2/c                -         Al1La1O3            #15       metal-oxide; Al1La1O3, ICSD #180415
A2B2C_mC20_15_ac_f_e-001.{ABC}                   mC20      C2/c                -         Cu2O2Pb1            #15       metal-oxide; Cu2O2Pb1, ICSD #400657
ABC4_mC24_15_a_e_2f-001.{ABC}                    mC24      C2/c                -         Hg1Mo1O4            #15       metal-oxide; Hg1Mo1O4, ICSD #2533
ABC4_mC24_15_e_e_2f-001.{ABC}                    mC24      C2/c                -         La1Nb1O4            #15       metal-oxide; La1Nb1O4, ICSD #10123
AB4C_mC24_15_e_2f_e-001.{ABC}                    mC24      C2/c                -         Nb1O4Y1             #15       metal-oxide; Nb1O4Y1, ICSD #20335
AB4C_mC24_15_a_2f_e-001.{ABC}                    mC24      C2/c                -         Hg1O4W1             #15       metal-oxide; Hg1O4W1, ICSD #169668
A2B4C_mC28_15_f_2f_e-002.{ABC}                   mC28      C2/c                -         Li2O4W1             #15       metal-oxide; Li2O4W1, ICSD #1044
A2B4C_mC28_15_f_2f_e-003.{ABC}                   mC28      C2/c                -         Li2O4W1             #15       metal-oxide; Li2O4W1, ICSD #14196
A2BC4_mC28_15_f_e_2f-001.{ABC}                   mC28      C2/c                -         Hg2Mo1O4            #15       metal-oxide; Hg2Mo1O4, ICSD #90084
A3B4C_mC32_15_cf_2f_e-001.{ABC}                  mC32      C2/c                -         Ag3O4V1             #15       metal-oxide; Ag3O4V1, ICSD #249417
A2B5C_mC32_15_f_e2f_c-001.{ABC}                  mC32      C2/c                -         Fe2O5Ti1            #15       metal-oxide; Fe2O5Ti1, ICSD #24416
ABC2_mC32_15_2e_2e_2f-001.{ABC}                  mC32      C2/c                -         Er1Na1O2            #15       metal-oxide; Er1Na1O2, ICSD #2739
A3B4C_oP16_31_ab_2ab_a-001.{ABC}                 oP16      Pmn2_{1}            -         Li3O4V1             #31       metal-oxide; Li3O4V1, ICSD #19002
ABC2_oP16_33_a_a_2a-002.{ABC}                    oP16      Pna2_{1}            -         Ga1Na1O2            #33       metal-oxide; Ga1Na1O2, ICSD #4416
ABC2_oP16_33_a_a_2a-003.ABC/BCA/CAB              oP16      Pna2_{1}            -         Ga1Li1O2            #33       metal-oxide; Ga1Li1O2, ICSD #18152, equivalent to GaLiO2 (http://aflow.org/prototype-encyclopedia/ABC2_oP16_33_a_a_2a.html, part 3)
ABC2_oP16_33_a_a_2a-004.{ABC}                    oP16      Pna2_{1}            -         Fe1Na1O2            #33       metal-oxide; Fe1Na1O2, ICSD #27117
ABC2_oP16_33_a_a_2a-005.{ABC}                    oP16      Pna2_{1}            -         Ga1Na1O2            #33       metal-oxide; Ga1Na1O2, ICSD #36652
ABC2_oP16_33_a_a_2a-006.{ABC}                    oP16      Pna2_{1}            -         Ag1Al1O2            #33       metal-oxide; Ag1Al1O2, ICSD #160643
AB8C2_oC22_35_a_ab3e_e-002.{ABC}                 oC22      Cmm2                -         Mo1O8V2             #35       metal-oxide; Mo1O8V2, ICSD #25378
ABC3_oC20_36_a_a_3a-001.{ABC}                    oC20      Cmc2_{1}            -         Ca1Ir1O3            #36       metal-oxide; Ca1Ir1O3, ICSD #159030
AB3C_oC20_36_a_ab_a-001.ABC/BCA/CAB              oC20      Cmc2_{1}            -         Mg1O3V1             #36       metal-oxide; Mg1O3V1, ICSD #284
A2BC2_oC20_36_2a_a_2a-001.{ABC}                  oC20      Cmc2_{1}            -         Na2Ni1O2            #36       metal-oxide; Na2Ni1O2, ICSD #14159
A2B3C_oC24_36_b_ab_a-003.{ABC}                   oC24      Cmc2_{1}            -         Cs2O3Pr1            #36       metal-oxide; Cs2O3Pr1, ICSD #1181
AB4C_oC24_36_a_4a_a-001.{ABC}                    oC24      Cmc2_{1}            -         La1O4Ta1            #36       metal-oxide; La1O4Ta1, ICSD #97688
AB5C2_oC32_36_a_5a_2a-001.{ABC}                  oC32      Cmc2_{1}            -         Li1O5V2             #36       metal-oxide; Li1O5V2, ICSD #23817
ABC3_oC10_38_a_b_ae-001.{ABC}                    oC10      Amm2                -         K1Nb1O3             #38       metal-oxide; K1Nb1O3, ICSD #9533
AB3C_oC10_38_b_ad_a-001.{ABC}                    oC10      Amm2                -         Ba1O3Ti1            #38       metal-oxide; Ba1O3Ti1, ICSD #154346
AB3C_oC10_38_b_ad_a-002.{ABC}                    oC10      Amm2                -         Ba1O3Ti1            #38       metal-oxide; Ba1O3Ti1, ICSD #161341
AB3C_oC10_38_a_ae_b-001.{ABC}                    oC10      Amm2                -         Ba1O3Ti1            #38       metal-oxide; Ba1O3Ti1, ICSD #161419
ABC_oC12_40_b_a_b-001.{ABC}                      oC12      Ama2                -         Cs1Cu1O1            #40       metal-oxide; Cs1Cu1O1, ICSD #37079
A2B4C_oC28_40_c_abc_b-001.{ABC}                  oC28      Ama2                -         Al2O4Pb1            #40       metal-oxide; Al2O4Pb1, ICSD #80128
A3BC_oI20_46_a2b_a_b-001.{ABC}                   oI20      Ima2                -         O3Sr1Ti1            #46       metal-oxide; O3Sr1Ti1, ICSD #182248
ABC2_oP4_47_a_g_s-001.{ABC}                      oP4       Pmmm                -         Au1K1O2             #47       metal-oxide; Au1K1O2, ICSD #15115
A3BC_oP5_47_bce_h_a-001.{ABC}                    oP5       Pmmm                -         O3Pb1Ti1            #47       metal-oxide; O3Pb1Ti1, ICSD #27949
AB4C_oP12_50_a_m_c-001.{ABC}                     oP12      Pban                -         Ce1O4Zr1            #50       metal-oxide; Ce1O4Zr1, ICSD #164740
A4BC2_oP14_55_gh_a_h-001.{ABC}                   oP14      Pbam                -         O4Pb1Sr2            #55       metal-oxide; O4Pb1Sr2, ICSD #4418
A2B4C_oP14_55_h_gh_a-001.{ABC}                   oP14      Pbam                -         Cd2O4Sn1            #55       metal-oxide; Cd2O4Sn1, ICSD #9010
AB2C_oP8_58_c_g_a-001.ABC/BCA/CAB                oP8       Pnnm                -         Li1O2Ru1            #58       metal-oxide; Li1O2Ru1, ICSD #48007
ABC2_oP8_59_a_a_2b-002.ABC/BCA/CAB               oP8       Pmmn                -         Mn1Na1O2            #59       metal-oxide; Mn1Na1O2, ICSD #16271
ABC2_oP8_59_a_a_2b-003.ABC/BCA/CAB               oP8       Pmmn                -         Li1Mn1O2            #59       metal-oxide; Li1Mn1O2, ICSD #84642
A2BC3_oP12_59_e_a_be-001.{ABC}                   oP12      Pmmn                -         Cu2Mg1O3            #59       metal-oxide; Cu2Mg1O3, ICSD #4202
A2B3C_oP12_59_e_be_a-001.{ABC}                   oP12      Pmmn                -         Cu2O3Sr1            #59       metal-oxide; Cu2O3Sr1, ICSD #95957
AB5C2_oP16_59_b_a2e_e-001.{ABC}                  oP16      Pmmn                -         Na1O5V2             #59       metal-oxide; Na1O5V2, ICSD #59345
AB2C_oP16_59_e_2e_e-001.{ABC}                    oP16      Pmmn                -         Li1O2Tb1            #59       metal-oxide; Li1O2Tb1, ICSD #21013
AB5C2_oP16_59_a_b2e_e-001.{ABC}                  oP16      Pmmn                -         Li1O5V2             #59       metal-oxide; Li1O5V2, ICSD #88643
ABC2_oP16_59_e_e_2e-001.{ABC}                    oP16      Pmmn                -         Gd1Li1O2            #59       metal-oxide; Gd1Li1O2, ICSD #27769
ABC_oP12_62_c_c_c-015.{ABC}                      oP12      Pnma                -         Li1O1Rb1            #62       metal-oxide; Li1O1Rb1, ICSD #65170
ABC2_oP16_62_c_c_2c-001.{ABC}                    oP16      Pnma                -         Eu1Li1O2            #62       metal-oxide; Eu1Li1O2, ICSD #422560
ABC_oC12_63_c_a_c-001.{ABC}                      oC12      Cmcm                -         Cs1Cu1O1            #63       metal-oxide; Cs1Cu1O1, ICSD #40161
ABC2_oC16_63_c_c_f-001.{ABC}                     oC16      Cmcm                -         Cu1K1O2             #63       metal-oxide; Cu1K1O2, ICSD #15095
ABC2_oC16_63_c_c_f-002.{ABC}                     oC16      Cmcm                -         Cs1Cu1O2            #63       metal-oxide; Cs1Cu1O2, ICSD #15097
ABC2_oC16_63_c_c_f-003.{ABC}                     oC16      Cmcm                -         Au1Na1O2            #63       metal-oxide; Au1Na1O2, ICSD #409547
ABC2_oC16_63_c_a_g-001.{ABC}                     oC16      Cmcm                -         Ba1Ni1O2            #63       metal-oxide; Ba1Ni1O2, ICSD #15760
ABC3_oC20_63_c_a_cf-002.{ABC}                    oC20      Cmcm                -         Ca1Ir1O3            #63       metal-oxide; Ca1Ir1O3, ICSD #159026
ABC3_oC20_63_a_c_cf-001.{ABC}                    oC20      Cmcm                -         Ir1Na1O3            #63       metal-oxide; Ir1Na1O3, ICSD #261371
AB4C_oC24_63_c_acf_c-001.{ABC}                   oC24      Cmcm                -         Bi1O4Re1            #63       metal-oxide; Bi1O4Re1, ICSD #10481
AB4C_oC24_63_a_fg_c-002.{ABC}                    oC24      Cmcm                -         In1O4V1             #63       metal-oxide; In1O4V1, ICSD #10431
A2B3C_oC24_63_g_ce_c-001.{ABC}                   oC24      Cmcm                -         Cs2O3Pb1            #63       metal-oxide; Cs2O3Pb1, ICSD #62140
AB4C_oC24_63_a_fg_c-003.{ABC}                    oC24      Cmcm                -         Cr1O4V1             #63       metal-oxide; Cr1O4V1, ICSD #27508
ABC4_oC24_63_c_a_fg-001.{ABC}                    oC24      Cmcm                -         Cr1Cu1O4            #63       metal-oxide; Cr1Cu1O4, ICSD #60825
ABC4_oC24_63_c_a_fg-002.{ABC}                    oC24      Cmcm                -         Cr1Hg1O4            #63       metal-oxide; Cr1Hg1O4, ICSD #416147
AB4C2_oC28_63_c_acf_f-001.{ABC}                  oC28      Cmcm                -         Ca1O4Tl2            #63       metal-oxide; Ca1O4Tl2, ICSD #79371
A2BC4_oC28_63_f_c_acf-001.{ABC}                  oC28      Cmcm                -         Al2Mg1O4            #63       metal-oxide; Al2Mg1O4, ICSD #161057
AB5C2_oC32_63_c_c2f_f-002.{ABC}                  oC32      Cmcm                -         Mg1O5V2             #63       metal-oxide; Mg1O5V2, ICSD #50979
AB5C2_oC32_63_c_c2f_f-003.{ABC}                  oC32      Cmcm                -         Li1O5V2             #63       metal-oxide; Li1O5V2, ICSD #50980
ABC_oC24_64_d_f_f-001.{ABC}                      oC24      Cmce                -         K1Li1O1             #64       metal-oxide; K1Li1O1, ICSD #30964
AB2C4_oC28_64_a_f_ef-001.{ABC}                   oC28      Cmce                -         Cu1La2O4            #64       metal-oxide; Cu1La2O4, ICSD #67836
AB2C4_oC28_64_a_f_ef-002.{ABC}                   oC28      Cmce                -         Cu1Nd2O4            #64       metal-oxide; Cu1Nd2O4, ICSD #86754
AB4C2_oC28_64_a_ef_f-001.{ABC}                   oC28      Cmce                -         Cu1O4Pr2            #64       metal-oxide; Cu1O4Pr2, ICSD #91073
ABC3_oC10_65_a_c_bf-001.{ABC}                    oC10      Cmmm                -         Na1Nb1O3            #65       metal-oxide; Na1Nb1O3, ICSD #28564
AB4C_oC12_65_a_hi_c-001.ABC/BCA/CAB              oC12      Cmmm                -         Co1O4Re1            #65       metal-oxide; Co1O4Re1, ICSD #72872
A2B3C_oC12_65_i_ai_c-001.{ABC}                   oC12      Cmmm                -         Cu2O3Sr1            #65       metal-oxide; Cu2O3Sr1, ICSD #99041
A2B3C_oC12_65_j_bj_a-001.{ABC}                   oC12      Cmmm                -         Li2O3Pr1            #65       metal-oxide; Li2O3Pr1, ICSD #154704 equivalent to Li2O3Pr (http://aflow.org/prototype-encyclopedia/A2B3C_oC12_65_h_bh_a.html, part 3)
AB3C4_oC16_65_b_ai_p-001.{ABC}                   oC16      Cmmm                -         Ba1Cu3O4            #65       metal-oxide; Ba1Cu3O4, ICSD #65881
AB2C_oC16_65_i_n_ad-001.{ABC}                    oC16      Cmmm                -         O1Ti2Zr1            #65       metal-oxide; O1Ti2Zr1, ICSD #9389
AB3C4_oC16_65_c_ai_p-001.{ABC}                   oC16      Cmmm                -         Ba1Cu3O4            #65       metal-oxide; Ba1Cu3O4, ICSD #83079
AB3C4_oC16_65_d_ai_p-001.{ABC}                   oC16      Cmmm                -         Ba1Cu3O4            #65       metal-oxide; Ba1Cu3O4, ICSD #89232
AB6C3_oC20_65_b_jp_af-001.{ABC}                  oC20      Cmmm                -         Mn1O6Pt3            #65       metal-oxide; Mn1O6Pt3, ICSD #35337
A3B5C2_oC20_65_ai_b2i_j-001.{ABC}                oC20      Cmmm                -         Cu6O10Sr4           #65       metal-oxide; Cu6O10Sr4, ICSD #50089
AB3C_oC20_65_i_acm_j-001.ABC/BCA/CAB             oC20      Cmmm                -         Mg1O3V1             #65       metal-oxide; Mg1O3V1, ICSD #15927
A3B7C_oC22_65_af_cjp_b-001.{ABC}                 oC22      Cmmm                -         La3O7Ta1            #65       metal-oxide; La3O7Ta1, ICSD #168909
A2BC4_oC28_65_ij_ac_ijo-001.{ABC}                oC28      Cmmm                -         La2Ni1O4            #65       metal-oxide; La2Ni1O4, ICSD #78338
A2BC_oC16_67_m_g_a-001.{ABC}                     oC16      Cmme                -         Ba2Na1O1            #67       metal-oxide; Ba2Na1O1, ICSD #411905
A2B3C2_oF28_69_h_ah_g-001.{ABC}                  oF28      Fmmm                -         Cu2O3Sr2            #69       metal-oxide; Cu2O3Sr2, ICSD #68676
AB4C2_oF28_69_a_cg_g-001.{ABC}                   oF28      Fmmm                -         Ni1O4Pr2            #69       metal-oxide; Ni1O4Pr2, ICSD #81577
A4BC_oF48_70_h_b_a-001.{ABC}                     oF48      Fddd                -         O4Tb1V1             #70       metal-oxide; O4Tb1V1, ICSD #88369
A2BC4_oF56_70_d_a_h-001.{ABC}                    oF56      Fddd                -         Cr2Mg1O4            #70       metal-oxide; Cr2Mg1O4, ICSD #290589
AB2C4_oF56_70_a_g_h-001.{ABC}                    oF56      Fddd                -         Mo1Na2O4            #70       metal-oxide; Mo1Na2O4, ICSD #151971
A2B2C_oI10_71_f_h_a-003.{ABC}                    oI10      Immm                -         K2O2Pd1             #71       metal-oxide; K2O2Pd1, ICSD #6158
A2B2C_oI10_71_e_f_a-001.{ABC}                    oI10      Immm                -         Li2O2Pd1            #71       metal-oxide; Li2O2Pd1, ICSD #19007
A2B2C_oI10_71_e_f_a-002.{ABC}                    oI10      Immm                -         Na2O2Pt1            #71       metal-oxide; Na2O2Pt1, ICSD #25018
A2B2C_oI10_71_e_f_a-003.{ABC}                    oI10      Immm                -         Ag2O2Pd1            #71       metal-oxide; Ag2O2Pd1, ICSD #51498
A2B2C_oI10_71_e_f_a-004.{ABC}                    oI10      Immm                -         Li2O2Pd1            #71       metal-oxide; Li2O2Pd1, ICSD #61199
AB2C2_oI10_71_a_e_f-001.{ABC}                    oI10      Immm                -         Cu1Li2O2            #71       metal-oxide; Cu1Li2O2, ICSD #25001
AB2C3_oI12_71_b_e_af-001.{ABC}                   oI12      Immm                -         Ba1Mn2O3            #71       metal-oxide; Ba1Mn2O3, ICSD #10038
A2BC3_oI12_71_e_a_de-002.{ABC}                   oI12      Immm                -         Ba2Cu1O3            #71       metal-oxide; Ba2Cu1O3, ICSD #68217
A2BC3_oI12_71_e_a_de-003.{ABC}                   oI12      Immm                -         Ca2Cu1O3            #71       metal-oxide; Ca2Cu1O3, ICSD #93651
ABC2_oI16_71_e_e_2f-001.{ABC}                    oI16      Immm                -         Li1Mn1O2            #71       metal-oxide; Li1Mn1O2, ICSD #15768
A2B4C3_oI18_71_f_n_ae-001.{ABC}                  oI18      Immm                -         Na2O4Pd3            #71       metal-oxide; Na2O4Pd3, ICSD #6157
A3B5C2_oI20_71_ae_def_e-001.{ABC}                oI20      Immm                -         Cu3O5Sr2            #71       metal-oxide; Cu3O5Sr2, ICSD #416904
A5B2C7_oI28_71_aef_g_dhn-001.{ABC}               oI28      Immm                -         Au5O2Rb7            #71       metal-oxide; Au5O2Rb7, ICSD #91309
A2B2C_oI20_72_j_j_a-002.{ABC}                    oI20      Ibam                -         K2O2Zn1             #72       metal-oxide; K2O2Zn1, ICSD #34603
AB3C2_oI24_72_c_bf_j-001.{ABC}                   oI24      Ibam                -         Ag1Na3O2            #72       metal-oxide; Ag1Na3O2, ICSD #16919
A3BC2_oI24_72_ce_b_j-001.{ABC}                   oI24      Ibam                -         Ag3Li1O2            #72       metal-oxide; Ag3Li1O2, ICSD #4204
ABC2_oI32_72_j_f_2j-001.{ABC}                    oI32      Ibam                -         Bi1Li1O2            #72       metal-oxide; Bi1Li1O2, ICSD #25385
A2BC_oI16_74_f_e_c-001.{ABC}                     oI16      Imma                -         O2Pb1Pd1            #74       metal-oxide; O2Pb1Pd1, ICSD #2277
AB3C_oI20_74_e_ef_c-001.{ABC}                    oI20      Imma                -         Ba1O3Pb1            #74       metal-oxide; Ba1O3Pb1, ICSD #15933
A2B2C_oI20_74_ac_i_e-001.{ABC}                   oI20      Imma                -         Cu2O2Tl1            #74       metal-oxide; Cu2O2Tl1, ICSD #36535
ABC4_oI24_74_e_e_hi-001.{ABC}                    oI24      Imma                -         Cr1Dy1O4            #74       metal-oxide; Cr1Dy1O4, ICSD #93789
A2BC4_tI28_79_c_2a_2c-001.{ABC}                  tI28      I4                  -         Bi2Cu1O4            #79       metal-oxide; Bi2Cu1O4, ICSD #12104
ABC_tI24_82_g_g_g-001.{ABC}                      tI24      I-4                 -         Ag1K1O1             #82       metal-oxide; Ag1K1O1, ICSD #25744
AB5C_tP14_85_a_cg_c-001.{ABC}                    tP14      P4/n                -         Mo1O5V1             #85       metal-oxide; Mo1O5V1, ICSD #27315
AB4C4_tI18_87_a_h_h-001.{ABC}                    tI18      I4/m                -         Ir1Na4O4            #87       metal-oxide; Ir1Na4O4, ICSD #67826
ABC_tI24_87_h_h_h-001.{ABC}                      tI24      I4/m                -         Cu1Na1O1            #87       metal-oxide; Cu1Na1O1, ICSD #15099
ABC_tI24_87_h_h_h-002.{ABC}                      tI24      I4/m                -         Cu1O1Rb1            #87       metal-oxide; Cu1O1Rb1, ICSD #15100
AB8C4_tI26_87_a_2h_h-001.{ABC}                   tI26      I4/m                -         K1O8Ru4             #87       metal-oxide; K1O8Ru4, ICSD #1562
AB8C4_tI26_87_a_2h_h-002.{ABC}                   tI26      I4/m                -         K2O16V8             #87       metal-oxide; K2O16V8, ICSD #100596
A2B5C_tI32_87_h_di_e-001.{ABC}                   tI32      I4/m                -         La4O10Re2           #87       metal-oxide; La4O10Re2, ICSD #81
A4BC_tI24_88_f_a_b-001.{ABC}                     tI24      I4_{1}/a            -         O4Pb1W1             #88       metal-oxide; O4Pb1W1, ICSD #16189
AB4C_tI24_88_a_f_b-001.{ABC}                     tI24      I4_{1}/a            -         Mo1O4Pb1            #88       metal-oxide; Mo1O4Pb1, ICSD #39137
AB4C_tI24_88_a_f_b-002.{ABC}                     tI24      I4_{1}/a            -         Ba1O4W1             #88       metal-oxide; Ba1O4W1, ICSD #52384
ABC4_tI24_88_a_b_f-003.{ABC}                     tI24      I4_{1}/a            -         Ca1Mo1O4            #88       metal-oxide; Ca1Mo1O4, ICSD #77334
AB4C_tI24_88_b_f_a-001.ABC/BCA/CAB               tI24      I4_{1}/a            -         Ba1O4W1             #88       metal-oxide; Ba1O4W1, ICSD #169094
A4B2C_tI28_88_f_d_a-001.{ABC}                    tI28      I4_{1}/a            -         O4Pd2Pr1            #88       metal-oxide; O4Pd2Pr1, ICSD #78843
A2BC4_tI28_88_d_a_f-001.{ABC}                    tI28      I4_{1}/a            -         Au2Ba1O4            #88       metal-oxide; Au2Ba1O4, ICSD #80327
ABC2_tP16_92_a_a_b-001.{ABC}                     tP16      P4_{1}2_{1}2        -         Al1Li1O2            #92       metal-oxide; Al1Li1O2, ICSD #23815
ABC2_tP16_92_a_a_b-002.{ABC}                     tP16      P4_{1}2_{1}2        -         Fe1Na1O2            #92       metal-oxide; Fe1Na1O2, ICSD #33768
ABC2_tP16_92_a_a_b-003.{ABC}                     tP16      P4_{1}2_{1}2        -         Al1K1O2             #92       metal-oxide; Al1K1O2, ICSD #151883
A3BC_tP5_99_ac_b_a-001.{ABC}                     tP5       P4mm                -         O3Pb1V1             #99       metal-oxide; O3Pb1V1, ICSD #152276
ABC3_tP5_99_b_a_ac-001.{ABC}                     tP5       P4mm                -         Bi1Fe1O3            #99       metal-oxide; Bi1Fe1O3, ICSD #188467
A6B5C_tI24_107_bc_ad_a-001.{ABC}                 tI24      I4mm                -         Na6O5Pb1            #107      metal-oxide; Na6O5Pb1, ICSD #15102
A2B4C_tI28_108_c_d_a-001.{ABC}                   tI28      I4cm                -         Bi2O4Pd1            #108      metal-oxide; Bi2O4Pd1, ICSD #9622
ABC2_tP8_113_c_a_e-001.{ABC}                     tP8       P-42_{1}m           -         Ca1Fe1O2            #113      metal-oxide; Ca1Fe1O2, ICSD #246244
ABC_tI24_119_g_i_i-001.{ABC}                     tI24      I-4m2               -         Cu1Li1O1            #119      metal-oxide; Cu1Li1O1, ICSD #282
ABC_tI24_119_g_i_i-002.{ABC}                     tI24      I-4m2               -         Ag1K1O1             #119      metal-oxide; Ag1K1O1, ICSD #24818
ABC_tI24_119_g_i_i-003.{ABC}                     tI24      I-4m2               -         Ag1Na1O1            #119      metal-oxide; Ag1Na1O1, ICSD #49752
ABC_tI24_119_g_i_i-004.{ABC}                     tI24      I-4m2               -         Ag1Cs1O1            #119      metal-oxide; Ag1Cs1O1, ICSD #49754
ABC_tI24_119_g_i_i-005.{ABC}                     tI24      I-4m2               -         Cu1Li1O1            #119      metal-oxide; Cu1Li1O1, ICSD #49755
A3B4C_tI16_121_bd_i_a-001.{ABC}                  tI16      I-42m               -         K3O4V1              #121      metal-oxide; K3O4V1, ICSD #4138
A3B4C_tI16_121_bd_i_a-002.{ABC}                  tI16      I-42m               -         Ag3O4V1             #121      metal-oxide; Ag3O4V1, ICSD #417470
A2BC6_tI18_121_e_a_di-001.{ABC}                  tI18      I-42m               -         La2Mo1O6            #121      metal-oxide; La2Mo1O6, ICSD #25611
AB3C8_tI24_121_a_bd_2i-001.{ABC}                 tI24      I-42m               -         Cr1K3O8             #121      metal-oxide; Cr1K3O8, ICSD #9356, equivalent to Cr1K3O8 (http://aflow.org/prototype-encyclopedia/AB3C8_tI24_121_a_bd_2i.html, part 3)
A3BC8_tI24_121_bd_a_2i-001.{ABC}                 tI24      I-42m               -         K3Nb1O8             #121      metal-oxide; K3Nb1O8, ICSD #30405
ABC2_tI16_122_b_a_d-001.{ABC}                    tI16      I-42d               -         Co1K1O2             #122      metal-oxide; Co1K1O2, ICSD #4199
A2BC4_tI28_122_d_a_e-001.{ABC}                   tI28      I-42d               -         Cr2Cu1O4            #122      metal-oxide; Cr2Cu1O4, ICSD #16708
AB3C_tP5_123_c_ae_b-001.{ABC}                    tP5       P4/mmm              -         Ba1O3Ti1            #123      metal-oxide; Ba1O3Ti1, ICSD #27965
AB3C_tP5_123_a_be_d-001.{ABC}                    tP5       P4/mmm              -         Ba1O3Ti1            #123      metal-oxide; Ba1O3Ti1, ICSD #109327
ABC3_tP5_123_c_b_ae-001.{ABC}                    tP5       P4/mmm              -         Ca1Mn1O3            #123      metal-oxide; Ca1Mn1O3, ICSD #168904
AB4C6_tP11_123_a_eh_bci-001.{ABC}                tP11      P4/mmm              -         Ba1Nb4O6            #123      metal-oxide; Ba1Nb4O6, ICSD #42006
A5B9C2_tP16_123_bfg_cegi_h-001.{ABC}             tP16      P4/mmm              -         Nb5O9Sr2            #123      metal-oxide; Nb5O9Sr2, ICSD #42007
AB3C_tP10_127_c_bg_a-001.{ABC}                   tP10      P4/mbm              -         Na1O3Ta1            #127      metal-oxide; Na1O3Ta1, ICSD #23322
ABC_tP6_129_c_a_c-002.{ABC}                      tP6       P4/nmm              -         K1Na1O1             #129      metal-oxide; K1Na1O1, ICSD #32743
AB4C2_tP14_131_c_n_j-001.{ABC}                   tP14      P4_{2}/mmc          -         Ca1O4Pt2            #131      metal-oxide; Ca1O4Pt2, ICSD #6159
A3BC2_tP12_136_ad_b_f-001.{ABC}                  tP12      P4_{2}/mnm          -         K3Ni1O2             #136      metal-oxide; K3Ni1O2, ICSD #262579
AB2C2_tI10_139_a_e_d-001.{ABC}                   tI10      I4/mmm              -         Bi1Ce2O2            #139      metal-oxide; Bi1Ce2O2, ICSD #9099 equivalent to Co2S2Tl (http://aflow.org/prototype-encyclopedia/A2B2C_tI10_139_d_e_a.TlCo2S2.html, part 3)
AB2C2_tI10_139_a_e_e-001.{ABC}                   tI10      I4/mmm              -         Hg1Na2O2            #139      metal-oxide; Hg1Na2O2, ICSD #25511
AB2C4_tI14_139_a_e_cd-001.{ABC}                  tI14      I4/mmm              -         Cu1Nd2O4            #139      metal-oxide; Cu1Nd2O4, ICSD #4203
AB2C4_tI14_139_a_e_cd-002.{ABC}                  tI14      I4/mmm              -         Cu1La2O4            #139      metal-oxide; Cu1La2O4, ICSD #261660
A2B4C_tI14_139_e_ce_a-001.{ABC}                  tI14      I4/mmm              -         Ba2O4Pb1            #139      metal-oxide; Ba2O4Pb1, ICSD #27113
AB2C4_tI14_139_a_e_ce-003.{ABC}                  tI14      I4/mmm              -         Cu1Nd2O4            #139      metal-oxide; Cu1Nd2O4, ICSD #86753
A3B2C6_tI22_139_ae_e_dg-001.{ABC}                tI22      I4/mmm              -         La3Ni2O6            #139      metal-oxide; La3Ni2O6, ICSD #249209
ABC_tI24_139_h_j_i-001.{ABC}                     tI24      I4/mmm              -         Ag1Na1O1            #139      metal-oxide; Ag1Na1O1, ICSD #40153
ABC_tI24_139_h_j_i-002.{ABC}                     tI24      I4/mmm              -         Ag1K1O1             #139      metal-oxide; Ag1K1O1, ICSD #40154
ABC_tI24_139_h_j_i-003.{ABC}                     tI24      I4/mmm              -         Cu1Li1O1            #139      metal-oxide; Cu1Li1O1, ICSD #40156
ABC_tI24_139_h_j_i-004.{ABC}                     tI24      I4/mmm              -         Ag1Cs1O1            #139      metal-oxide; Ag1Cs1O1, ICSD #40160
A7B3C2_tI24_139_aeg_be_e-001.{ABC}               tI24      I4/mmm              -         O7Sr3Ti2            #139      metal-oxide; O7Sr3Ti2, ICSD #20294, equivalent to O7Sr3Ti2 (http://aflow.org/prototype-encyclopedia/A7B3C2_tI24_139_aeg_be_e.html, part 3)
ABC_tI24_139_h_i_j-001.{ABC}                     tI24      I4/mmm              -         Cu1O1Rb1            #139      metal-oxide; Cu1O1Rb1, ICSD #40159
ABC_tI24_139_h_i_j-002.{ABC}                     tI24      I4/mmm              -         Au1O1Rb1            #139      metal-oxide; Au1O1Rb1, ICSD #409552
A4B3C8_tI30_139_2e_ae_cdg-001.{ABC}              tI30      I4/mmm              -         Nd4Ni3O8            #139      metal-oxide; Nd4Ni3O8, ICSD #51097
A3BC_tI20_140_ah_b_c-001.{ABC}                   tI20      I4/mcm              -         O3Sr1Zr1            #140      metal-oxide; O3Sr1Zr1, ICSD #1522
AB3C_tI20_140_c_ah_b-001.{ABC}                   tI20      I4/mcm              -         Hf1O3Sr1            #140      metal-oxide; Hf1O3Sr1, ICSD #164621
ABC2_tI16_141_b_a_e-001.{ABC}                    tI16      I4_{1}/amd          -         Gd1Na1O2            #141      metal-oxide; Gd1Na1O2, ICSD #22040
ABC2_tI16_141_a_b_e-002.ABC/BCA/CAB              tI16      I4_{1}/amd          -         Fe2Li2O4            #141      metal-oxide; Fe2Li2O4, ICSD #31149
AB2C2_tI20_141_a_d_e-001.{ABC}                   tI20      I4_{1}/amd          -         Ba1Cu2O2            #141      metal-oxide; Ba1Cu2O2, ICSD #9456
AB4C_tI24_141_a_h_b-001.{ABC}                    tI24      I4_{1}/amd          -         Cs1O4Re1            #141      metal-oxide; Cs1O4Re1, ICSD #72817
A2B2C3_tI28_141_d_c_ae-001.{ABC}                 tI28      I4_{1}/amd          -         Ag2Cu2O3            #141      metal-oxide; Ag2Cu2O3, ICSD #51502
A2B4C_tI28_141_d_h_a-001.{ABC}                   tI28      I4_{1}/amd          -         Mn2O4Zn1            #141      metal-oxide; Mn2O4Zn1, ICSD #15305
AB2C4_tI28_141_a_d_h-001.{ABC}                   tI28      I4_{1}/amd          -         Cd1Mn2O4            #141      metal-oxide; Cd1Mn2O4, ICSD #24258
A2BC4_tI28_141_d_a_h-001.{ABC}                   tI28      I4_{1}/amd          -         Cr2Ni1O4            #141      metal-oxide; Cr2Ni1O4, ICSD #37023
A2BC4_tI28_141_d_a_h-002.{ABC}                   tI28      I4_{1}/amd          -         Cr2Fe1O4            #141      metal-oxide; Cr2Fe1O4, ICSD #183968
A2BC4_tI28_141_d_a_h-003.{ABC}                   tI28      I4_{1}/amd          -         Cr2Cu1O4            #141      metal-oxide; Cr2Cu1O4, ICSD #246080
ABC2_tI32_141_c_d_h-001.ABC/BCA/CAB              tI32      I4_{1}/amd          -         Li1Mn1O2            #141      metal-oxide; Li1Mn1O2, ICSD #40486
AB2C_tI32_141_d_h_c-001.ABC/BCA/CAB              tI32      I4_{1}/amd          -         Li1O2Ti1            #141      metal-oxide; Li1O2Ti1, ICSD #164158
ABC3_hR5_146_a_a_b-001.{ABC}                     hR5       R3                  -         Bi1Fe1O3            #146      metal-oxide; Bi1Fe1O3, ICSD #162264
A3B9C4_hR16_146_3a_3b_4a-001.{ABC}               hR16      R3                  -         Ba3O9Yb4            #146      metal-oxide; Ba3O9Yb4, ICSD #33239
A3B9C4_hR16_146_3a_3b_4a-002.{ABC}               hR16      R3                  -         Ba3O9Y4             #146      metal-oxide; Ba3O9Y4, ICSD #87118
A3B4C9_hR16_146_3a_4a_3b-001.{ABC}               hR16      R3                  -         Ba3Ho4O9            #146      metal-oxide; Ba3Ho4O9, ICSD #33807
AB8C2_hP11_147_a_dg_d-001.{ABC}                  hP11      P-3                 -         Mn1O8Re2            #147      metal-oxide; Mn1O8Re2, ICSD #51014
AB8C2_hP11_147_a_dg_d-002.{ABC}                  hP11      P-3                 -         Ni1O8Re2            #147      metal-oxide; Ni1O8Re2, ICSD #51016
AB3C_hR10_148_c_f_c-002.ABC/BCA/CAB              hR10      R-3                 -         Cd1O3Ti1            #148      metal-oxide; Cd1O3Ti1, ICSD #15989
AB3C_hR10_148_c_f_c-003.{ABC}                    hR10      R-3                 -         Fe1O3Ti1            #148      metal-oxide; Fe1O3Ti1, ICSD #153674
A8B6C_hR15_148_cf_f_a-001.{ABC}                  hR15      R-3                 -         Li8O6Sn1            #148      metal-oxide; Li8O6Sn1, ICSD #1180, equivalent to Li7O6Sn1 (partially occupied Li, http://aflow.org/prototype-encyclopedia/A8B6C_hR15_148_cf_f_a.html, part 3)
AB2C_hP12_152_b_c_a-001.{ABC}                    hP12      P3_{1}21            -         Ba1O2Zn1            #152      metal-oxide; Ba1O2Zn1, ICSD #25559
AB2C_hP12_152_a_c_b-001.{ABC}                    hP12      P3_{1}21            -         Ba1O2Zn1            #152      metal-oxide; Ba1O2Zn1, ICSD #25812
AB2C_hP12_154_a_c_b-001.{ABC}                    hP12      P3_{2}21            -         Hg1O2Sr1            #154      metal-oxide; Hg1O2Sr1, ICSD #69739
AB3C_hR5_160_a_b_a-001.{ABC}                     hR5       R3m                 -         Ba1O3Ti1            #160      metal-oxide; Ba1O3Ti1, ICSD #6102
A2B3C2_hR7_160_2a_b_2a-001.{ABC}                 hR7       R3m                 -         K2O3Sn2             #160      metal-oxide; K2O3Sn2, ICSD #2216
ABC3_hR10_161_a_a_b-002.{ABC}                    hR10      R3c                 -         Na1Nb1O3            #161      metal-oxide; Na1Nb1O3, ICSD #9645
ABC3_hR10_161_a_a_b-003.{ABC}                    hR10      R3c                 -         Li1Nb1O3            #161      metal-oxide; Li1Nb1O3, ICSD #28296
ABC3_hR10_161_a_a_b-004.{ABC}                    hR10      R3c                 -         Ga1La1O3            #161      metal-oxide; Ga1La1O3, ICSD #51036
AB3C_hR10_161_a_b_a-001.{ABC}                    hR10      R3c                 -         Ni1O3Pb1            #161      metal-oxide; Ni1O3Pb1, ICSD #187684
A2B3C_hR12_161_2a_b_a-001.{ABC}                  hR12      R3c                 -         Li2O3Re1            #161      metal-oxide; Li2O3Re1, ICSD #200999
A5B6C2_hP13_162_ef_k_d-001.{ABC}                 hP13      P-31m               -         Ag5O6Pb2            #162      metal-oxide; Ag5O6Pb2, ICSD #40058
A4B7C2_hP13_162_de_ak_c-001.{ABC}                hP13      P-31m               -         In4O7Rb2            #162      metal-oxide; In4O7Rb2, ICSD #6321
AB3C2_hP6_164_b_e_d-001.{ABC}                    hP6       P-3m1               -         Cs1Cu3O2            #164      metal-oxide; Cs1Cu3O2, ICSD #413342
A2B8C_hP11_164_d_di_a-001.{ABC}                  hP11      P-3m1               -         Mo2O8Zr1            #164      metal-oxide; Mo2O8Zr1, ICSD #59999
A7B4C2_hP13_164_ai_2d_c-001.{ABC}                hP13      P-3m1               -         O7Tl4V2             #164      metal-oxide; O7Tl4V2, ICSD #72809
A3B8C3_hP14_164_ad_di_bd-001.{ABC}               hP14      P-3m1               -         K3O8V3              #164      metal-oxide; K3O8V3, ICSD #100782
A4B2C_hP14_164_di_abd_d-001.{ABC}                hP14      P-3m1               -         O4Tl2W1             #164      metal-oxide; O4Tl2W1, ICSD #8212
ABC2_hR4_166_a_b_c-005.{ABC}                     hR4       R-3m                -         Ag1Cr1O2            #166      metal-oxide; Ag1Cr1O2, ICSD #4149
ABC2_hR4_166_a_b_c-006.{ABC}                     hR4       R-3m                -         Cu1La1O2            #166      metal-oxide; Cu1La1O2, ICSD #18102
ABC2_hR4_166_a_b_c-007.{ABC}                     hR4       R-3m                -         Cu1Eu1O2            #166      metal-oxide; Cu1Eu1O2, ICSD #18106
ABC2_hR4_166_a_b_c-008.{ABC}                     hR4       R-3m                -         Cu1In1O2            #166      metal-oxide; Cu1In1O2, ICSD #55687
ABC2_hR4_166_b_a_c-001.{ABC}                     hR4       R-3m                -         Al1Cu1O2            #166      metal-oxide; Al1Cu1O2, ICSD #25593
A2BC_hR4_166_c_a_b-001.{ABC}                     hR4       R-3m                -         O2Rb1Sc1            #166      metal-oxide; O2Rb1Sc1, ICSD #31959
AB2C_hR4_166_a_c_b-003.{ABC}                     hR4       R-3m                -         Hg1O2Sr1            #166      metal-oxide; Hg1O2Sr1, ICSD #165068
ABC2_hR4_166_b_a_c-002.ABC/BCA/CAB               hR4       R-3m                -         Mo1Na1O2            #166      metal-oxide; Mo1Na1O2, ICSD #166514
ABC3_hR5_166_b_a_d-001.{ABC}                     hR5       R-3m                -         Bi1Fe1O3            #166      metal-oxide; Bi1Fe1O3, ICSD #20288
A2BC2_hR5_166_c_a_c-001.{ABC}                    hR5       R-3m                -         Ag2Ni1O2            #166      metal-oxide; Ag2Ni1O2, ICSD #160574
A2B4C_hR7_166_c_2c_a-001.{ABC}                   hR7       R-3m                -         Fe2O4Yb1            #166      metal-oxide; Fe2O4Yb1, ICSD #4192
A2B3C2_hR7_166_ab_d_c-001.{ABC}                  hR7       R-3m                -         K2O3Sn2             #166      metal-oxide; K2O3Sn2, ICSD #15511, equivalent to Ni3Pb2S2 (shandite, http://aflow.org/prototype-encyclopedia/A3B2C2_hR7_166_d_ab_c.html, part 3)
A2B4C_hR7_166_c_2c_a-002.{ABC}                   hR7       R-3m                -         Fe2O4Y1             #166      metal-oxide; Fe2O4Y1, ICSD #67701
A2BC4_hR7_166_c_a_2c-001.{ABC}                   hR7       R-3m                -         Fe2In1O4            #166      metal-oxide; Fe2In1O4, ICSD #157323
AB2C_hR8_166_c_ad_c-001.{ABC}                    hR8       R-3m                -         Al1O2Tl1            #166      metal-oxide; Al1O2Tl1, ICSD #29010
AB3C2_hR12_166_c_h_ae-001.{ABC}                  hR12      R-3m                -         K1O3Pd2             #166      metal-oxide; K1O3Pd2, ICSD #248051
AB4C8_hR13_166_a_bd_ch-001.{ABC}                 hR13      R-3m                -         Ba1Ni4O8            #166      metal-oxide; Ba1Ni4O8, ICSD #20898
AB4C8_hR13_166_a_bd_ch-002.{ABC}                 hR13      R-3m                -         Ca1Ni4O8            #166      metal-oxide; Ca1Ni4O8, ICSD #40470
A3B2C8_hR13_166_ac_c_ch-001.{ABC}                hR13      R-3m                -         Ba3Cr2O8            #166      metal-oxide; Ba3Cr2O8, ICSD #9457
A2B8C3_hR13_166_c_ch_e-001.{ABC}                 hR13      R-3m                -         Ca2O8Pt3            #166      metal-oxide; Ca2O8Pt3, ICSD #65412
A4B8C_hR13_166_bd_ch_a-001.{ABC}                 hR13      R-3m                -         Ni4O8Sr1            #166      metal-oxide; Ni4O8Sr1, ICSD #40469
A8B3C2_hR13_166_ch_ac_c-001.{ABC}                hR13      R-3m                -         O8Pb3V2             #166      metal-oxide; O8Pb3V2, ICSD #27651
AB4C2_hR14_166_c_ch_ad-001.{ABC}                 hR14      R-3m                -         Al1O4V2             #166      metal-oxide; Al1O4V2, ICSD #151457
AB3C_hR15_166_bc_dh_ac-001.{ABC}                 hR15      R-3m                -         Ba1O3Ru1            #166      metal-oxide; Ba1O3Ru1, ICSD #10253
A3BC4_hR16_167_e_a_be-001.{ABC}                  hR16      R-3c                -         Ba3Ni1O4            #167      metal-oxide; Ba3Ni1O4, ICSD #30662
AB2C4_hP14_182_b_f_cg-001.{ABC}                  hP14      P6_{3}22            -         Ba1Fe2O4            #182      metal-oxide; Ba1Fe2O4, ICSD #2769
ABC3_hP10_186_b_a_c-001.{ABC}                    hP10      P6_{3}mc            -         Ba1Ni1O3            #186      metal-oxide; Ba1Ni1O3, ICSD #15761
AB3C_hP10_186_a_a2b_b-001.{ABC}                  hP10      P6_{3}mc            -         La1O3Tl1            #186      metal-oxide; La1O3Tl1, ICSD #200088
AB2C4_hP14_186_b_ab_bc-001.{ABC}                 hP14      P6_{3}mc            -         Mn1Na2O4            #186      metal-oxide; Mn1Na2O4, ICSD #39504
AB2C_hP8_187_bc_gh_i-001.{ABC}                   hP8       P-6m2               -         Cs1O2Y1             #187      metal-oxide; Cs1O2Y1, ICSD #49652
ABC3_hP10_187_ad_i_jk-001.{ABC}                  hP10      P-6m2               -         Ba1Co1O3            #187      metal-oxide; Ba1Co1O3, ICSD #88670
A3B4C_hP16_190_af_cg_d-001.{ABC}                 hP16      P-62c               -         Be3O4Sr1            #190      metal-oxide; Be3O4Sr1, ICSD #26179
A2B3C_hP6_191_d_f_a-001.{ABC}                    hP6       P6/mmm              -         Mn2O3Ta1            #191      metal-oxide; Mn2O3Ta1, ICSD #15995
A2BC_hP8_194_f_c_a-001.{ABC}                     hP8       P6_{3}/mmc          -         O2Rb1Sc1            #194      metal-oxide; O2Rb1Sc1, ICSD #1270
ABC2_hP8_194_c_a_f-002.{ABC}                     hP8       P6_{3}/mmc          -         Ag1Fe1O2            #194      metal-oxide; Ag1Fe1O2, ICSD #2786, equivalent to hexagonal delafossite (AlCuO2, http://aflow.org/prototype-encyclopedia/ABC2_hP8_194_a_c_f.html, part 3)
ABC2_hP8_194_c_a_f-003.{ABC}                     hP8       P6_{3}/mmc          -         Cs1Nd1O2            #194      metal-oxide; Cs1Nd1O2, ICSD #27336
ABC2_hP8_194_c_a_f-004.{ABC}                     hP8       P6_{3}/mmc          -         Cu1Fe1O2            #194      metal-oxide; Cu1Fe1O2, ICSD #66546
ABC2_hP8_194_a_c_f-002.{ABC}                     hP8       P6_{3}/mmc          -         Na1Nb1O2            #194      metal-oxide; Na1Nb1O2, ICSD #29282
ABC2_hP8_194_a_c_f-003.{ABC}                     hP8       P6_{3}/mmc          -         Al1Au1O2            #194      metal-oxide; Al1Au1O2, ICSD #95663
AB2C_hP8_194_c_f_a-001.{ABC}                     hP8       P6_{3}/mmc          -         Cu1O2Y1             #194      metal-oxide; Cu1O2Y1, ICSD #35580
AB2C_hP8_194_c_f_a-002.{ABC}                     hP8       P6_{3}/mmc          -         Cu1O2Sc1            #194      metal-oxide; Cu1O2Sc1, ICSD #60847
AB3C_hP10_194_c_bf_a-001.{ABC}                   hP10      P6_{3}/mmc          -         Al1O3Y1             #194      metal-oxide; Al1O3Y1, ICSD #27100
ABC3_hP10_194_c_a_bf-001.{ABC}                   hP10      P6_{3}/mmc          -         Ga1In1O3            #194      metal-oxide; Ga1In1O3, ICSD #30339
ABC3_hP10_194_a_c_bf-001.{ABC}                   hP10      P6_{3}/mmc          -         In1Mn1O3            #194      metal-oxide; In1Mn1O3, ICSD #67671
AB3C_hP10_194_c_bf_a-002.{ABC}                   hP10      P6_{3}/mmc          -         Mn1O3Y1             #194      metal-oxide; Mn1O3Y1, ICSD #73361
ABC2_hP16_194_ac_f_ef-001.{ABC}                  hP16      P6_{3}/mmc          -         Ag1Co1O2            #194      metal-oxide; Ag1Co1O2, ICSD #261608
A2B3C2_cI28_199_a_b_a-001.{ABC}                  cI28      I2_{1}3             -         K2O3Pb2             #199      metal-oxide; K2O3Pb2, ICSD #1412
A4B2C7_cF52_216_e_bc_ag-001.{ABC}                cF52      F-43m               -         Al4Cu2O7            #216      metal-oxide; Al4Cu2O7, ICSD #100355
A4B3C_cP16_223_e_c_a.{ABC}                       cP16      Pm-3n               -         O4Pd3Sr1            #223      metal-oxide; O4Pd3Sr1, ICSD #16537
ABC3_cF40_225_c_ab_e-001.{ABC}                   cF40      Fm-3m               -         Ba2Bi2O6            #225      metal-oxide; Ba2Bi2O6, ICSD #68714
A3B6C_cF40_225_bc_e_a-001.{ABC}                  cF40      Fm-3m               -         Eu3O6Ta1            #225      metal-oxide; Eu3O6Ta1, ICSD #4148
A3B6C_cF40_225_bc_e_a-002.{ABC}                  cF40      Fm-3m               -         Ba3O6W1             #225      metal-oxide; Ba3O6W1, ICSD #76437
A6BC8_cF60_225_d_a_ce-001.{ABC}                  cF60      Fm-3m               -         Mg6Mn1O8            #225      metal-oxide; Mg6Mn1O8, ICSD #24710
A6B8C_cF60_225_d_f_a-001.{ABC}                   cF60      Fm-3m               -         Cu6O8Pb1            #225      metal-oxide; Cu6O8Pb1, ICSD #280596
A4B3C_cF64_225_f_d_ab-001.{ABC}                  cF64      Fm-3m               -         O4Pd3Tl1            #225      metal-oxide; O4Pd3Tl1, ICSD #2275
AB2C_cF64_227_c_e_d-001.ABC/BCA/CAB              cF64      Fd-3m               -         Li1O2Ti1            #227      metal-oxide; Li1O2Ti1, ICSD #48128

*****************************************************************************************************************************************************
ANRL - metal carbides - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2B_aP12_2_4i_2i-001.{AB}                        aP12      P-1                 -         C2Ca1               #2        metal-carbide; C2Ca1, ICSD #66663
AB4_mP10_11_e_4e-001.{AB}                        mP10      P2_{1}/m            -         C1Fe4               #11       metal-carbide; C1Fe4, ICSD #187143 equivalent to 1271.AB,1271.BA [HTQC]
A2B_mC12_12_2i_i-003.{AB}                        mC12      C2/m                -         C2Ca1               #12       metal-carbide; C2Ca1, ICSD #54185 equivalent to 1244.AB,1244.BA [HTQC]
A2B_mC12_12_2i_i-004.{AB}                        mC12      C2/m                -         C2Ca1               #12       metal-carbide; C2Ca1, ICSD #54188 equivalent to 1245.AB,1245.BA [HTQC], equivalent to CaC2-III (http://aflow.org/prototype-encyclopedia/A2B_mC12_12_2i_i.html, part 3)
A2B_mC12_12_2i_i-005.{AB}                        mC12      C2/m                -         C2Ir1               #12       metal-carbide; C2Ir1, ICSD #181488 equivalent to 1246.AB,1246.BA [HTQC]
A5B6_mC22_12_agh_ij-001.{AB}                     mC22      C2/m                -         C5Nb6               #12       metal-carbide; C5Nb6, ICSD #20695 equivalent to 1279.AB,1279.BA [HTQC]
A4B_mP10_14_2e_a-001.{AB}                        mP10      P2_{1}/c            -         C4Ir1               #14       metal-carbide; C4Ir1, ICSD #181493 equivalent to 1276.AB,1276.BA [HTQC]
A2B_mP12_14_2e_e-010.{AB}                        mP12      P2_{1}/c            -         C2Ca1               #14       metal-carbide; C2Ca1, ICSD #24074 equivalent to 1289.AB,1289.BA [HTQC]
A2B_mC12_15_f_e-001.{AB}                         mC12      C2/c                -         C2Ca1               #15       metal-carbide; C2Ca1, ICSD #54184 equivalent to 1208.AB,1208.BA [HTQC]
A3B5_oP16_55_ah_d2g-001.{AB}                     oP16      Pbam                -         C3Ir5               #55       metal-carbide; C3Ir5, ICSD #181485 equivalent to 1317.AB,1317.BA [HTQC]
A3B2_oP10_58_ag_g-001.{AB}                       oP10      Pnnm                -         C3Mg2               #58       metal-carbide; C3Mg2, ICSD #71941 equivalent to 1275.AB,1275.BA [HTQC]
AB2_oP12_62_c_2c-002.{AB}                        oP12      Pnma                -         C1Nb2               #62       metal-carbide; C1Nb2, ICSD #31973 equivalent to 1285.AB,1285.BA [HTQC]
AB2_oP12_62_c_2c-003.{AB}                        oP12      Pnma                -         C1Fe2               #62       metal-carbide; C1Fe2, ICSD #187138 equivalent to 1282.AB,1282.BA [HTQC]
AB_oP16_62_2c_2c-001.{AB}                        oP16      Pnma                -         C2Rb2               #62       metal-carbide; C2Rb2, ICSD #51529 equivalent to 1215.AB,1215.BA,1313.AB,1313.BA [HTQC]
A3B_oP16_62_3c_c-001.{AB}                        oP16      Pnma                -         C3Ir1               #62       metal-carbide; C3Ir1, ICSD #181492 equivalent to 1316.AB,1316.BA [HTQC]
AB3_oP16_62_c_cd-003.{AB}                        oP16      Pnma                -         C1Fe3               #62       metal-carbide; C1Fe3, ICSD #167667 equivalent to 1307.AB,1307.BA [HTQC]
AB3_oP16_62_c_cd-004.{AB}                        oP16      Pnma                -         C1Cr3               #62       metal-carbide; C1Cr3, ICSD #617486
A4B_oC20_68_i_a-001.{AB}                         oC20      Ccce                -         C4Ir1               #68       metal-carbide; C4Ir1, ICSD #181496 equivalent to 1277.AB,1277.BA [HTQC]
AB_oI8_71_g_e-002.{AB}                           oI8       Immm                -         C2Li2               #71       metal-carbide; C2Li2, ICSD #25705 equivalent to 1220.AB,1220.BA [HTQC]
AB4_tI10_87_a_h-002.{AB}                         tI10      I4/m                -         C1Fe4               #87       metal-carbide; C1Fe4, ICSD #187144 equivalent to 1236.AB,1236.BA [HTQC]
A2B_tP3_123_h_a-001.{AB}                         tP3       P4/mmm              -         C2Re1               #123      metal-carbide; C2Re1, ICSD #184664 equivalent to 1229.AB,1229.BA [HTQC]
A2B_tP6_136_f_a-003.{AB}                         tP6       P4_{2}/mnm          -         C2Mg1               #136      metal-carbide; C2Mg1, ICSD #88057 equivalent to 1248.AB,1248.BA [HTQC]
A2B_tI12_140_h_a-002.{AB}                        tI12      I4/mcm              -         C2Ir1               #140      metal-carbide; C2Ir1, ICSD #181489 equivalent to 1247.AB,1247.BA [HTQC]
AB_tI32_142_f_e-001.{AB}                         tI32      I4_{1}/acd          -         C2Na2               #142      metal-carbide; C2Na2, ICSD #28066 equivalent to 1315.AB,1315.BA [HTQC]
AB_tI32_142_f_e-002.{AB}                         tI32      I4_{1}/acd          -         C2K2                #142      metal-carbide; C2K2, ICSD #36142 equivalent to 1314.AB,1314.BA [HTQC]
A4B_hP15_152_2ac_b-001.{AB}                      hP15      P3_{1}21            -         C4Ir1               #152      metal-carbide; C4Ir1, ICSD #181495
A4B3_hR7_160_4a_3a-001.{AB}                      hR7       R3m                 -         Al4C3               #160      metal-carbide; Al4C3, ICSD #14397 equivalent to 1254.AB,1254.BA [HTQC]
A5B8_hR13_160_2ab_2a2b-001.{AB}                  hR13      R3m                 -         C5Ti8               #160      metal-carbide; C5Ti8, ICSD #618924
AB2_hP3_164_a_d-002.{AB}                         hP3       P-3m1               -         C1Sc2               #164      metal-carbide; C1Sc2, ICSD #280743 equivalent to 1227.AB,1227.BA [HTQC]
AB2_hR3_166_a_c-002.{AB}                         hR3       R-3m                -         C1Y2                #166      metal-carbide; C1Y2, ICSD #22283 equivalent to 1207.AB,1207.BA [HTQC]
AB2_hR3_166_a_c-003.{AB}                         hR3       R-3m                -         Ba1C2               #166      metal-carbide; Ba1C2, ICSD #186576 equivalent to 1225.AB,1225.BA [HTQC]
AB_hR4_166_c_c-004.{AB}                          hR4       R-3m                -         C1Ru1               #166      metal-carbide; C1Ru1, ICSD #188285 equivalent to 1233.AB,1233.BA [HTQC]
A5B8_hR13_166_abd_ch-001.{AB}                    hR13      R-3m                -         C5Ti8               #166      metal-carbide; C5Ti8, ICSD #20822 equivalent to 1294.AB,1294.BA [HTQC]
AB_hP4_186_b_b-005.AB                            hP4       P6_{3}mc            -         C1Sc1               #186      metal-carbide; C1Sc1, ICSD #189087 equivalent to 1230.AB [HTQC]
AB2_hP3_187_a_i-001.{AB}                         hP3       P-6m2               -         C1Mo2               #187      metal-carbide; C1Mo2, ICSD #43669
AB2_hP3_191_a_c-001.{AB}                         hP3       P6/mmm              -         C1Fe2               #191      metal-carbide; C1Fe2, ICSD #162103 equivalent to 1241.AB,1241.BA [HTQC]
A8B_hP9_191_cl_b-001.{AB}                        hP9       P6/mmm              -         C8Cs1               #191      metal-carbide; C8Cs1, ICSD #74641 equivalent to 1270.AB,1270.BA [HTQC]
A2B_hP6_194_e_c-002.{AB}                         hP6       P6_{3}/mmc          -         C2Os1               #194      metal-carbide; C2Os1, ICSD #168280 equivalent to 1221.AB,1221.BA [HTQC]
A2B_hP6_194_f_c-002.{AB}                         hP6       P6_{3}/mmc          -         C2Re1               #194      metal-carbide; C2Re1, ICSD #184660 equivalent to 1249.AB,1249.BA [HTQC]
A6B_hP14_194_i_c-001.{AB}                        hP14      P6_{3}/mmc          -         C6Eu1               #194      metal-carbide; C6Eu1, ICSD #169041 equivalent to 1302.AB,1302.BA [HTQC]
A6B_hP14_194_i_c-002.{AB}                        hP14      P6_{3}/mmc          -         C6Yb1               #194      metal-carbide; C6Yb1, ICSD #601565 equivalent to 1303.AB,1303.BA [HTQC]
AB_cP8_198_a_a-003.{AB}                          cP8       P2_{1}3             -         C1Os1               #198      metal-carbide; C1Os1, ICSD #168277 equivalent to 1266.AB,1266.BA [HTQC]
A3B4_cI28_220_a_c-001.{AB}                       cI28      I-43d               -         C3Sc4               #220      metal-carbide; C3Sc4, ICSD #42760 equivalent to 1301.AB,1301.BA [HTQC]
A2B_cP12_224_e_b-001.{AB}                        cP12      Pn-3m               -         C2Ca1               #224      metal-carbide; C2Ca1, ICSD #31092 equivalent to 1290.AB,1290.BA [HTQC]
AB_cF64_227_cd_e-001.{AB}                        cF64      Fd-3m               -         C1Ti1               #227      metal-carbide; C1Ti1, ICSD #618950 equivalent to 1312.AB,1312.BA [HTQC]

*****************************************************************************************************************************************************
ANRL - metal carbides - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2BC_mC16_9_2a_a_a-001.{ABC}                     mC16      Cc                  -         C2Co1Nd1            #9        metal-carbide; C2Co1Nd1, ICSD #55570
A2BC_mC16_9_2a_a_a-002.{ABC}                     mC16      Cc                  -         C2Co1Nd1            #9        metal-carbide; C2Co1Nd1, ICSD #67375
A2BC_mC16_9_2a_a_a-003.{ABC}                     mC16      Cc                  -         C2Co1Pr1            #9        metal-carbide; C2Co1Pr1, ICSD #658525
A3B2C2_mC14_12_ai_i_i-002.{ABC}                  mC14      C2/m                -         C3Cr2Ho2            #12       metal-carbide; C3Cr2Ho2, ICSD #62083
A3B2C2_mC14_12_ai_i_i-003.{ABC}                  mC14      C2/m                -         C3Ho2Mo2            #12       metal-carbide; C3Ho2Mo2, ICSD #88511
A3B2C2_mC14_12_ai_i_i-004.{ABC}                  mC14      C2/m                -         C3Ce2Mo2            #12       metal-carbide; C3Ce2Mo2, ICSD #417827
A3B2C2_mC14_12_ai_i_i-005.{ABC}                  mC14      C2/m                -         C3Er2Mo2            #12       metal-carbide; C3Er2Mo2, ICSD #617669
A5B4C3_mC24_12_b2i_2i_ai-001.{ABC}               mC24      C2/m                -         C5Ca4Ni3            #12       metal-carbide; C5Ca4Ni3, ICSD #71440
A4BC3_mC32_12_2j_i_ghi-001.{ABC}                 mC32      C2/m                -         C4Co1Sc3            #12       metal-carbide; C4Co1Sc3, ICSD #236391
A4BC3_mC32_12_2j_i_bc2i-001.{ABC}                mC32      C2/m                -         C4Ru1Sc3            #12       metal-carbide; C4Ru1Sc3, ICSD #420074
A5B2C4_oP11_25_c2g_g_adh-001.{ABC}               oP11      Pmm2                -         C5Ni2Yb4            #25       metal-carbide; C5Ni2Yb4, ICSD #73156
A2BC_oC8_38_d_b_a-001.{ABC}                      oC8       Amm2                -         C2Nd1Rh1            #38       metal-carbide; C2Nd1Rh1, ICSD #63064
A2BC_oP16_57_2d_d_c-001.{ABC}                    oP16      Pbcm                -         C2Cs1Na1            #57       metal-carbide; C2Cs1Na1, ICSD #189824
A2BC_oP16_59_2e_ab_e-001.{ABC}                   oP16      Pmmn                -         C2Cr1Sc1            #59       metal-carbide; C2Cr1Sc1, ICSD #80372
A2BC_oP16_62_2c_c_c-002.{ABC}                    oP16      Pnma                -         C2Dy1Mo1            #62       metal-carbide; C2Dy1Mo1, ICSD #617606
A2BC_oP16_62_d_c_c-002.{ABC}                     oP16      Pnma                -         C2Cs1K1             #62       metal-carbide; C2Cs1K1, ICSD #189822
A2BC_oP16_62_d_c_c-003.{ABC}                     oP16      Pnma                -         C2Cs1Rb1            #62       metal-carbide; C2Cs1Rb1, ICSD #189823
A2BC_oC16_63_g_c_a-001.{ABC}                     oC16      Cmcm                -         C2Gd1Ru1            #63       metal-carbide; C2Gd1Ru1, ICSD #80312
A2B2C_oC20_63_ac_f_c-001.{ABC}                   oC20      Cmcm                -         C2Cr2V1             #63       metal-carbide; C2Cr2V1, ICSD #20297
A4BC3_oC16_65_p_a_bi-001.{ABC}                   oC16      Cmmm                -         C4Os1Sc3            #65       metal-carbide; C4Os1Sc3, ICSD #420076
A7B2C5_oC28_65_aegj_i_bq-001.{ABC}               oC28      Cmmm                -         C7Re2Sc5            #65       metal-carbide; C7Re2Sc5, ICSD #72382
A4BC3_oI16_71_n_b_af-001.{ABC}                   oI16      Immm                -         C4Fe1Sc3            #71       metal-carbide; C4Fe1Sc3, ICSD #72863
A4B2C_oI28_72_2j_j_a-001.{ABC}                   oI28      Ibam                -         C4Er2Fe1            #72       metal-carbide; C4Er2Fe1, ICSD #42968
A2BC_tP16_76_2a_a_a-002.{ABC}                    tP16      P4_{1}              -         C2Ce1Rh1            #76       metal-carbide; C2Ce1Rh1, ICSD #67327
A2BC_tP4_123_g_a_d-001.{ABC}                     tP4       P4/mmm              -         C2Cu1Rb1            #123      metal-carbide; C2Cu1Rb1, ICSD #391118
AB2C_tP4_123_a_g_d-001.{ABC}                     tP4       P4/mmm              -         Ag1C2K1             #123      metal-carbide; Ag1C2K1, ICSD #410874
AB2C_tP4_123_a_g_d-002.{ABC}                     tP4       P4/mmm              -         Au1C2Na1            #123      metal-carbide; Au1C2Na1, ICSD #411254
AB2C2_tP5_123_b_e_ac-001.{ABC}                   tP5       P4/mmm              -         C1Co2Mn2            #123      metal-carbide; C1Co2Mn2, ICSD #44353
A2BC_tP8_129_2c_a_c-001.{ABC}                    tP8       P4/nmm              -         C2Co1Sc1            #129      metal-carbide; C2Co1Sc1, ICSD #62598
AB2C_tP8_131_a_j_f-001.{ABC}                     tP8       P4_{2}/mmc          -         Ag1C2Cs1            #131      metal-carbide; Ag1C2Cs1, ICSD #410873
A2BC_tP8_131_j_a_f-001.{ABC}                     tP8       P4_{2}/mmc          -         C2Cu1K1             #131      metal-carbide; C2Cu1K1, ICSD #412038
A4B5C2_hR11_160_4a_5a_2a-001.{ABC}               hR11      R3m                 -         Al4C5Zr2            #160      metal-carbide; Al4C5Zr2, ICSD #173676
A4B6C3_hR13_160_4a_6a_3a-001.{ABC}               hR13      R3m                 -         Al4C6Zr3            #160      metal-carbide; Al4C6Zr3, ICSD #173677
A2B2C_hP5_164_c_d_a-001.{ABC}                    hP5       P-3m1               -         C2Na2Pd1            #164      metal-carbide; C2Na2Pd1, ICSD #50172
A2BC2_hP5_164_c_a_d-001.{ABC}                    hP5       P-3m1               -         C2Pd1Rb2            #164      metal-carbide; C2Pd1Rb2, ICSD #94394
A4B5C2_hR11_166_2c_a2c_c-001.{ABC}               hR11      R-3m                -         Al4C5Hf2            #166      metal-carbide; Al4C5Hf2, ICSD #161586
A4B6C3_hR13_166_2c_3c_ac-001.{ABC}               hR13      R-3m                -         Al4C6Hf3            #166      metal-carbide; Al4C6Hf3, ICSD #161585
A3B3C_hP14_186_3b_3b_a-001.{ABC}                 hP14      P6_{3}mc            -         Al3C3Sc1            #186      metal-carbide; Al3C3Sc1, ICSD #43477
AB2C_hP4_187_a_g_d-001.{ABC}                     hP4       P-6m2               -         Ag1C2Li1            #187      metal-carbide; Ag1C2Li1, ICSD #410868
ABC2_hP8_194_a_c_f-004.{ABC}                     hP8       P6_{3}/mmc          -         C1Pb1Ti2            #194      metal-carbide; C1Pb1Ti2, ICSD #42926
A2BC_hP8_194_f_c_a-002.{ABC}                     hP8       P6_{3}/mmc          -         C2Cr1Sc1            #194      metal-carbide; C2Cr1Sc1, ICSD #80373
ABC2_hP8_194_c_a_f-005.{ABC}                     hP8       P6_{3}/mmc          -         Al1C1W2             #194      metal-carbide; Al1C1W2, ICSD #165101
ABC2_hP8_194_c_a_f-006.{ABC}                     hP8       P6_{3}/mmc          -         Al1C1Ta2            #194      metal-carbide; Al1C1Ta2, ICSD #181247
AB2C3_hP12_194_b_f_af-001.{ABC}                  hP12      P6_{3}/mmc          -         Al1C2Ti3            #194      metal-carbide; Al1C2Ti3, ICSD #93503
A2BC3_hP12_194_f_b_af-001.{ABC}                  hP12      P6_{3}/mmc          -         C2Sn1Ti3            #194      metal-carbide; C2Sn1Ti3, ICSD #160572
A3B3C_hP14_194_cf_df_a-001.{ABC}                 hP14      P6_{3}/mmc          -         Al3C3Sc1            #194      metal-carbide; Al3C3Sc1, ICSD #62308
AB3C4_hP16_194_c_ae_2f-002.{ABC}                 hP16      P6_{3}/mmc          -         Al1C3Ta4            #194      metal-carbide; Al1C3Ta4, ICSD #157843
AB3C4_hP16_194_c_ae_2f-003.{ABC}                 hP16      P6_{3}/mmc          -         Al1C3Ta4            #194      metal-carbide; Al1C3Ta4, ICSD #159456

*****************************************************************************************************************************************************
ANRL - carbo-nitrides - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A11B4_oP15_16_acgqtu_u-001.{AB}                  oP15      P222                -         C11N4               #16       carbo-nitride; C11N4, ICSD #184897
AB2_oC12_36_a_2a-001.{AB}                        oC12      Cmc2_{1}            -         C1N2                #36       carbo-nitride; C1N2, ICSD #247680
AB_oP16_61_c_c-002.{AB}                          oP16      Pbca                -         C2N2                #61       carbo-nitride; C2N2, ICSD #15870
A11B4_tP15_111_abcmn_n-001.{AB}                  tP15      P-42m               -         C11N4               #111      carbo-nitride; C11N4, ICSD #184896
AB2_tP6_113_a_e-001.{AB}                         tP6       P-42_{1}m           -         C1N2                #113      carbo-nitride; C1N2, ICSD #247676
AB2_tI6_119_a_f-001.{AB}                         tI6       I-4m2               -         C1N2                #119      carbo-nitride; C1N2, ICSD #247678
AB2_tI24_122_d_e-001.{AB}                        tI24      I-42d               -         C1N2                #122      carbo-nitride; C1N2, ICSD #247677
A3B4_hR7_160_b_ab-001.{AB}                       hR7       R3m                 -         C3N4                #160      carbo-nitride; C3N4, ICSD #41952
AB2_hP6_164_c_2d-001.{AB}                        hP6       P-3m1               -         C1N2                #164      carbo-nitride; C1N2, ICSD #247679
A3B4_hP14_176_h_ch-001.{AB}                      hP14      P6_{3}/m            -         C3N4                #176      carbo-nitride; C3N4, ICSD #41950, equivalent to beta-Si3N4 (http://aflow.org/prototype-encyclopedia/A4B3_hP14_176_ch_h.html, part 3)
A3B4_hP14_187_jk_adjk-001.{AB}                   hP14      P-6m2               -         C3N4                #187      carbo-nitride; C3N4, ICSD #83265
A3B4_cP7_215_c_e-001.{AB}                        cP7       P-43m               -         C3N4                #215      carbo-nitride; C3N4, ICSD #41951 equivalent to 707.AB,707.BA,1255.AB,1255.BA [HTQC]
A3B4_cI28_220_b_c-001.{AB}                       cI28      I-43d               -         C3N4                #220      carbo-nitride; C3N4, ICSD #83263

*****************************************************************************************************************************************************
ANRL - metal carbo-nitrides - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB2C2_aP15_2_ai_3i_3i-001.{ABC}                  aP15      P-1                 -         C1N2Tl2             #2        metal carbo-nitride; C1N2Tl2, ICSD #417297
A4B3C_aP16_2_4i_3i_i-001.{ABC}                   aP16      P-1                 -         C4N3Na1             #2        metal carbo-nitride; C4N3Na1, ICSD #31926
A4BC3_aP16_2_4i_i_3i-002.{ABC}                   aP16      P-1                 -         C4K1N3              #2        metal carbo-nitride; C4K1N3, ICSD #77046
ABC_mC6_8_a_a_a-001.ABC/CAB/BCA                  mC6       Cm                  -         C1K1N1              #8        metal carbo-nitride; C1K1N1, ICSD #27350
ABC_mC12_9_a_a_a-001.{ABC}                       mC12      Cc                  -         C1K1N1              #9        metal carbo-nitride; C1K1N1, ICSD #173942
AB2C2_mC10_12_a_i_i-002.{ABC}                    mC10      C2/m                -         C1K2N2              #12       metal carbo-nitride; C1K2N2, ICSD #411094
AB2C2_mC10_12_a_i_i-003.{ABC}                    mC10      C2/m                -         C1N2Na2             #12       metal carbo-nitride; C1N2Na2, ICSD #411341
ABC2_mP16_14_e_e_2e-002.{ABC}                    mP16      P2_{1}/c            -         C1Hg1N2             #14       metal carbo-nitride; C1Hg1N2, ICSD #412278
AB2C_oP16_33_a_2a_a-001.{ABC}                    oP16      Pna2_{1}            -         C1N2Pb1             #33       metal carbo-nitride; C1N2Pb1, ICSD #16600
AB4C3_oI32_46_b_2bc_bc-001.{ABC}                 oI32      Ima2                -         Ag1C4N3             #46       metal carbo-nitride; Ag1C4N3, ICSD #43823
ABC_oP6_59_a_a_b-001.ABC/CAB/BCA                 oP6       Pmmn                -         C1N1Na1             #59       metal carbo-nitride; C1N1Na1, ICSD #77172
ABC_oP12_62_c_c_c-016.{ABC}                      oP12      Pnma                -         C1Li1N1             #62       metal carbo-nitride; C1Li1N1, ICSD #77321
AB2C_oP16_62_c_2c_c-002.{ABC}                    oP16      Pnma                -         C1N2Sr1             #62       metal carbo-nitride; C1N2Sr1, ICSD #75040
AB2C_oP16_62_c_2c_c-003.{ABC}                    oP16      Pnma                -         C1N2Pb1             #62       metal carbo-nitride; C1N2Pb1, ICSD #410915
ABC2_oC16_63_c_a_f-001.{ABC}                     oC16      Cmcm                -         C1Cu1N2             #63       metal carbo-nitride; C1Cu1N2, ICSD #161460
ABC2_tI16_122_a_b_d-002.{ABC}                    tI16      I-42d               -         C1Mg1N2             #122      metal carbo-nitride; C1Mg1N2, ICSD #44110
AB2C_tI32_122_d_e_d-001.{ABC}                    tI32      I-42d               -         C1N2Zn1             #122      metal carbo-nitride; C1N2Zn1, ICSD #280523
AB2C2_tI10_139_a_d_e-002.{ABC}                   tI10      I4/mmm              -         C1Li2N2             #139      metal carbo-nitride; C1Li2N2, ICSD #200369, equivalent to CLi2N2 (http://aflow.org/prototype-encyclopedia/AB2C2_tI10_139_a_d_e.Li2CN2.html, part 3)
A3B2C6_hR11_155_d_c_f-001.{ABC}                  hR11      R32                 -         C3Lu2N6             #155      metal carbo-nitride; C3Lu2N6, ICSD #240311
ABC_hR3_160_a_a_a-002.ABC/BAC/CAB                hR3       R3m                 -         Ag1C1N1             #160      metal carbo-nitride; Ag1C1N1, ICSD #85783
ABC2_hR4_160_a_a_2a-001.{ABC}                    hR4       R3m                 -         C1Cd1N2             #160      metal carbo-nitride; C1Cd1N2, ICSD #95265
A6B3C2_hR11_160_6a_3a_2a-001.{ABC}               hR11      R3m                 -         Al6C3N2             #160      metal carbo-nitride; Al6C3N2, ICSD #14399
A8B3C4_hR15_160_8a_3a_4a-001.{ABC}               hR15      R3m                 -         Al8C3N4             #160      metal carbo-nitride; Al8C3N4, ICSD #14401
ABC2_hR4_166_b_a_c-003.{ABC}                     hR4       R-3m                -         C1Ca1N2             #166      metal carbo-nitride; C1Ca1N2, ICSD #25763
ABC2_hR4_166_b_a_c-004.{ABC}                     hR4       R-3m                -         C1Ca1N2             #166      metal carbo-nitride; C1Ca1N2, ICSD #31100
ABC2_hR4_166_b_a_c-005.{ABC}                     hR4       R-3m                -         C1Mg1N2             #166      metal carbo-nitride; C1Mg1N2, ICSD #75039
A6B3C2_hR11_166_3c_ac_c-001.{ABC}                hR11      R-3m                -         Al6C3N2             #166      metal carbo-nitride; Al6C3N2, ICSD #41260
AB2C_hR4_166_a_c_b-004.{ABC}                     hR4       R-3m                -         C1N2Sr1             #166      metal carbo-nitride; C1N2Sr1, ICSD #59860
A8B3C4_hR15_166_4c_ac_2c-001.{ABC}               hR15      R-3m                -         Al8C3N4             #166      metal carbo-nitride; Al8C3N4, ICSD #41261
AB2C_hP8_194_c_f_a-003.{ABC}                     hP8       P6_{3}/mmc          -         C1N2Ni1             #194      metal carbo-nitride; C1N2Ni1, ICSD #249388
A2BC2_cP10_215_e_ab_e-001.ABC/CAB/BCA            cP10      P-43m               -         C2Cd1N2             #215      metal carbo-nitride; C2Cd1N2, ICSD #20748

*************************************************************************************************************************************************************
ANRL PROTOTYPES - Part 3
*************************************************************************************************************************************************************

ANRL Label                                                           Pearson   SG symbol           Struk         Prototype               SG #      notes
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A5B11CD8E_aP26_1_5a_11a_a_8a_a-001.{ABCDE}                           aP26      P1                  -             C5H11NaO8S              #1        NaC5H11O8S (part 3) 
A2B2C5_aP18_2_2i_2i_5i-001.{ABC}                                     aP18      P-1                 -             B2Co2O5                 #2        Co2B2O5 (part 3) 
ABC8D3_aP26_2_i_i_8i_3i-001.{ABCD}                                   aP26      P-1                 S6_{8}        AlNaO8Si3               #2        Albite (NaAlSi3O8) (part 3) 
AB3C3_aP28_2_2i_6i_6i-001.{ABC}                                      aP28      P-1                 G5_{1}        BH3O3                   #2        Boric Acid (H3BO3) (part 3) 
A3B_aP32_2_12i_4i-001.{AB}                                           aP32      P-1                 -             O3W                     #2        delta-WO3 (part 3) 
A3B2C10D3_aP36_2_ah2i_2i_10i_3i-001.{ABCD}                           aP36      P-1                 -             Co3H2O10Se3             #2        Co3(SeO3)3.H2O (part 3) 
AB10C9D_aP42_2_ae_10i_9i_i-001.{ABCD}                                aP42      P-1                 H4_{10}       CuH10O9S                #2        Chalcanthite (CuSO4.5H2O) (part 3) 
A2B7C2_aP44_2_4i_14i_4i-001.{ABC}                                    aP44      P-1                 -             Ho2O7Si2                #2        alpha-Ho2Si2O7 (part 3) 
A12B2CD12_aP54_2_12i_2i_i_12i-001.{ABCD}                             aP54      P-1                 -             H12N2NiO12              #2        Ni(NO3)2(H2O)6 (part 3) 
A2B2C5D_mP20_4_2a_2a_5a_a-001.{ABCD}                                 mP20      P2_{1}              H4_{8}        H2Li2O5S                #4        Li2SO4.H2O (part 3) 
A3B6C_mP20_4_3a_6a_a-001.{ABC}                                       mP20      P2_{1}              -             Ca3O6U                  #4        Ca3UO6 (part 3) 
A11B2C2_mP60_4_22a_4a_4a-001.{ABC}                                   mP60      P2_{1}              -             O11P2W2                 #4        W2O3(PO4)2 (part 3) 
ABC3_mC10_5_b_a_ac-001.{ABC}                                         mC10      C2                  -             (Ba,Ca)CO3              #5        C2 (Ba,Ca)CO3 (part 3) 
A2B_mC12_5_2c_c-001.{AB}                                             mC12      C2                  -             As2Nb                   #5        NbAs2 (part 3) 
AB3_mC16_5_c_3c-001.{AB}                                             mC16      C2                  D0_{15}       AlCl3                   #5        D015 (AlCl3) (obsolete) (part 3) 
AB6C18D4E2_mC62_5_a_2b2c_9c_2c_c-001.{ABCDE}                         mC62      C2                  -             CaCu6O18P4Rb2           #5        Rb2CaCu6(PO4)4O2 (part 3) 
A2B2C9D2_mC90_5_ab2c_3c_b13c_3c-001.{ABCD}                           mC90      C2                  H4_{7}        Ca2H2O9S2               #5        Bassanite [CaSO4(H2O)0.5] (part 3) 
AB2_mP12_7_2a_4a-001.{AB}                                            mP12      Pc                  -             AuTe2                   #7        Calaverite (AuTe2) (part 3) 
A6B2C15D4_mP54_7_6a_2a_15a_4a-001.{ABCD}                             mP54      Pc                  -             Ca6Na2O15Si4            #7        Na2Ca6Si4O15 (part 3) 
A8B23_mP124_7_16a_46a-001.{AB}                                       mP124     Pc                  -             Mo8O23                  #7        Low-Temperature Mo8O23 (part 3) 
ABC2_mC8_8_a_a_b-001.{ABC}                                           mC8       Cm                  F5_{11}       KNO2                    #8        F511 (KNO2) (obsolete) (part 3) 
A13B4_mC102_8_17a11b_8a2b-001.{AB}                                   mC102     Cm                  -             Al13Co4                 #8        Monoclinic Co4Al13 (part 3) equivalent to 493,494 [HTQC]
A3B5C4D2_mC56_9_3a_5a_4a_2a-001.{ABCD}                               mC56      Cc                  -             Mg3O5(OH)4Si2           #9        Chrysotile (Mg3Si2O5(OH)4) (part 3) 
A2B4C9D2_mC68_9_2a_4a_9a_2a-001.{ABCD}                               mC68      Cc                  S5_{4}        Al2H4O9Si2              #9        Nacrite [Al2Si2O5(OH)4] (part 3) 
A6B36C11_mC212_9_6a_36a_11a-001.{ABC}                                mC212     Cc                  -             Cs6O36W11               #9        Cs6W11O36 (part 3) 
ABC_mP6_11_e_e_e-001.{ABC}                                           mP6       P2_{1}/m            -             HO2Y                    #11       O(OH)Y (part 3) 
A3B_mP8_11_3e_e-001.{AB}                                             mP8       P2_{1}/m            -             Se3Zr                   #11       ZrSe3 (part 3) 
A2B5C2_mP18_11_2e_e2f_2e-001.{ABC}                                   mP18      P2_{1}/m            K0_{1}        K2O5S2                  #11       K2S2O5 (part 3) 
AB2CD6_mP20_11_e_2e_e_2e2f-001.{ABCD}                                mP20      P2_{1}/m            -             BaC2CaO6                #11       Barytocalcite (BaCa(CO3)2) (part 3) 
A7B2C2_mP22_11_3e2f_2e_ab-001.{ABC}                                  mP22      P2_{1}/m            -             O7Si2Y2                 #11       y-Y2Si2O7 (part 3) 
A2B_mC6_12_i_a-001.{AB}                                              mC6       C2/m                -             Cl2Cu                   #12       Tolbachite (CuCl2) (part 3) 
A3B4_mC14_12_ai_2i-001.{AB}                                          mC14      C2/m                D7_{a}        Ni3Sn4                  #12       delta-Ni3Sn4 (part 3) 
AB2C_mC16_12_g_2i_i-001.{ABC}                                        mC16      C2/m                -             FeSe2Tl                 #12       Monoclinic FeTlSe2 (part 3) 
AB2_mC18_12_ai_3i-002.{AB}                                           mC18      C2/m                -             NbTe2                   #12       NbTe2 (part 3) 
AB6C2_mC18_12_a_3i_i-004.{ABC}                                       mC18      C2/m                -             PdSe6Ta2                #12       Ta2PdSe6 (part 3) 
ABC3_mC20_12_g_i_ij-001.{ABC}                                        mC20      C2/m                -             MnPS3                   #12       MnPS3 (part 3) 
A3BC2_mC24_12_ij_h_gi-001.{ABC}                                      mC24      C2/m                B36           H3LiO2                  #12       LiOH.H2O (part 3) 
AB_mC24_12_3i_3i-001.{AB}                                            mC24      C2/m                -             AsSi                    #12       SiAs (part 3) 
A13B4_mC34_12_b6i_2i-001.{AB}                                        mC34      C2/m                -             Al13Os4                 #12       Os4Al13 (part 3) 
AB6C2D_mC40_12_ad_gh4i_j_bc-001.ABCD/DACB/CABD/                      mC40      C2/m                -             NiO6Sr2Te               #12       Sr2NiTeO6 (part 3) 
DCBA/BCAD/DBAC/CDAB/ADBC/BADC/CBDA/ACDB/BDCA
A2B12CD6_mC42_12_i_2i2j_a_ij-001.{ABCD}                              mC42      C2/m                J1_{7}        Cl2H12MgO6              #12       Bischofite (MgCl2.6H2O) (part 3) 
AB5_mC48_12_2i_ac5i2j-001.{AB}                                       mC48      C2/m                D2_{2}        MgZn5                   #12       D22 (MgZn5?) (problematic) (part 3) 
AB8C4_mC52_12_i_gi3j_2j-001.{ABC}                                    mC52      C2/m                S6_{7}        AlKO8Si3                #12       Sanidine (KAlSi3O8) (part 3) 
A2B2C5D24E8_mC82_12_h_i_agh_2i5j_2j-001.{ABCDE}                      mC82      C2/m                S4_{2}        Ca2H2Mg5O24Si8          #12       Tremolite (Ca2Mg5Si8O22(OH)2) (part 3) 
A5B2C10D2E2_mC84_12_acghj_bdi_5j_2i_j-001.{ABCDE}                    mC84      C2/m                -             Al5Fe2O10(OH)2Si2       #12       Staurolite (Al5Fe2O10(OH)2Si2) (part 3) 
A13B4_mC102_12_dg8i5j_4ij-001.{AB}                                   mC102     C2/m                -             Al13Fe4                 #12       Al13Fe4 (part 3) 
A8B2CD15E2_mC112_12_2i3j_j_ad_g4i5j_2i-001.{ABCDE}                   mC112     C2/m                H4_{23}       H8K2MnO12S2             #12       Manganese-leonite [K2Mn(SO4)2.4H2O] (part 3) 
AB6C11D6E4_mC112_12_e_gi2j_i5j_2i2j_2j-001.{ABCDE}                   mC112     C2/m                S4_{5}        (H2O)(OH)6O11Mg6Si4     #12       Chrysotile (H4Mg3Si2O9) (part 3) 
A8B23_mP62_13_4g_c11g-001.{AB}                                       mP62      P2/c                -             Mo8O23                  #13       High-Temperature Mo8O23 (part 3) 
ABC_mP12_14_e_e_e-002.{ABC}                                          mP12      P2_{1}/c            E0_{7}        AsFeS                   #14       Arsenopyrite (FeAsS) (part 3) 
A2B3C2_mP14_14_e_ae_e-001.{ABC}                                      mP14      P2_{1}/c            -             Cl2Hg3O2                #14       HgCl2.2HgO (part 3) 
AB3_mP16_14_e_3e-002.{AB}                                            mP16      P2_{1}/c            -             beta-B2H6               #14       beta-B2H6 (part 3) 
ABC2_mP16_14_e_e_2e-003.{ABC}                                        mP16      P2_{1}/c            -             KNO2                    #14       KNO2 III (part 3) 
ABCD_mP16_14_e_e_e_e-001.{ABCD}                                      mP16      P2_{1}/c            -             ClCuHO                  #14       Cu(OH)Cl (part 3) 
AB_mP16_14_2e_2e-001.{AB}                                            mP16      P2_{1}/c            -             ClI                     #14       alpha-ICl (part 3) 
AB3_mP16_14_e_3e-003.{AB}                                            mP16      P2_{1}/c            -             B2H6                    #14       B2H6 (P21/c) (part 3) 
ABC2_mP16_14_e_e_2e-004.{ABC}                                        mP16      P2_{1}/c            E0_{6}        MnHO2                   #14       Manganite (gamma-MnO(OH)) (part 3) 
AB_mP16_14_2e_2e-002.{AB}                                            mP16      P2_{1}/c            -             AsLi                    #14       LiAs (part 3) 
A2B3_mP20_14_2e_3e-002.{AB}                                          mP20      P2_{1}/c            D5_{f}        As2S3                   #14       Orpiment (As2S3) (part 3) 
A3BC_mP20_14_3e_e_e-001.{ABC}                                        mP20      P2_{1}/c            -             Cl3CuK                  #14       Sanguite (KCuCl3) (part 3) 
AB6C2D_mP20_14_a_3e_e_d-001.ABCD/DACB/CABD/                          mP20      P2_{1}/c            -             MnO6Sr2Te               #14       Sr2MnTeO6 (part 3) 
DCBA/BCAD/DBAC/CDAB/ADBC/BADC/CBDA/ACDB/BDCA
AB6C3_mP20_14_a_3e_de-001.{ABC}                                      mP20      P2_{1}/c            J2_{6}        AlF6Na3                 #14       Cryolite (Na3AlF6) (part 3) 
A2B5C4_mP22_14_e_c2e_2e-001.{ABC}                                    mP22      P2_{1}/c            -             Cl2O5Sb4                #14       Sb4O5Cl2 (part 3) 
A4B2C4D_mP22_14_2e_e_2e_a-001.{ABCD}                                 mP22      P2_{1}/c            -             C4K2N4Ni                #14       K2Ni(CN)4 (part 3) 
A9B2_mP22_14_a4e_e-001.{AB}                                          mP22      P2_{1}/c            D8_{d}        Co2Al9                  #14       Co2Al9 (part 3) 
A4BC_mP24_14_4e_e_e-001.{ABC}                                        mP24      P2_{1}/c            -             O7Si2Y2                 #14       gamma-Y2Si2O7 (part 3) 
AB4C_mP24_14_ab_4e_e-001.{ABC}                                       mP24      P2_{1}/c            -             AuBr4K                  #14       Anhydrous KAuBr4 (part 3) 
ABCD3_mP24_14_e_e_e_3e-001.{ABCD}                                    mP24      P2_{1}/c            G0_{12}       CHNaO3                  #14       Nahcolite (NaHCO3) (part 3) 
AB_mP24_14_3e_3e-001.{AB}                                            mP24      P2_{1}/c            -             C6Cl6                   #14       epsilon-1,2,3,4,5,6-Hexachlorocyclohexane (C6Cl6) (part 3) 
AB4C_mP24_14_e_4e_e-001.{ABC}                                        mP24      P2_{1}/c            K4_{1}        (NH4)O4S                #14       Ammonium Persulfate [(NH4)2S2O8] (part 3) 
AB4C_mP24_14_e_4e_e-002.{ABC}                                        mP24      P2_{1}/c            -             LaO4P                   #14       Monazite (LaPO4) (part 3) 
ABC4_mP24_14_e_e_4e-001.{ABC}                                        mP24      P2_{1}/c            H0_{9}        AgMnO4                  #14       AgMnO4 (part 3) 
A2B4C_mP28_14_abe_4e_e-001.{ABC}                                     mP28      P2_{1}/c            -             Cu2O4Se                 #14       Monoclinic Cu2OSeO3 (part 3) 
A4BCD_mP28_14_4e_e_e_e-001.{ABCD}                                    mP28      P2_{1}/c            H0_{10}       Cl4IK.H2O               #14       KICl4.H2O (part 3) 
A2B3C2D8_mP30_14_e_ce_e_4e-001.{ABCD}                                mP30      P2_{1}/c            G7_{4}        C2Cu3H2O8               #14       Azurite [Cu3(CO3)2(OH)2] (part 3) 
A2B5C_mP32_14_2e_5e_ab-001.{ABC}                                     mP32      P2_{1}/c            -             Ca2O5U                  #14       Ca2UO5 (part 3) 
A2B5C_mP32_14_2e_5e_e-001.{ABC}                                      mP32      P2_{1}/c            -             Gd2O5Si                 #14       Gd2SiO5 (RE2SiO5 X1) (part 3) 
A3B_mP32_14_6e_2e-001.{AB}                                           mP32      P2_{1}/c            -             O3W                     #14       gamma-WO3 (part 3) 
AB4C2D_mP32_14_e_4e_2e_e-001.{ABCD}                                  mP32      P2_{1}/c            H4_{19}       AuBr4(H2O)2K            #14       KAuBr4.2H2O (part 3) 
AB_mP32_14_4e_4e-001.{AB}                                            mP32      P2_{1}/c            -             AsS                     #14       Pararealgar (AsS) (part 3) 
AB_mP32_14_4e_4e-002.{AB}                                            mP32      P2_{1}/c            B_{l}         AsS                     #14       Realgar (AsS) (part 3) 
A7B2C_mP40_14_7e_2e_e-001.{ABC}                                      mP40      P2_{1}/c            K6_{2}        F7K2Nb                  #14       K2NbF7 (part 3) 
A6B4C2D6E2FG6_mP54_14_3e_2e_e_3e_e_a_3e-001.{ABCDEFG}                mP54      P2_{1}/c            -             C6H4K2N6O2PtS6          #14       K2Pt(SCN)6.2H2O (part 3) 
A11B3_mP56_14_11e_3e-001.{AB}                                        mP56      P2_{1}/c            -             Cs11O3                  #14       Cs11O3 (part 3) 
AB3C_mP60_14_3e_9e_3e-001.{ABC}                                      mP60      P2_{1}/c            S3_{3} (II)   CaO3Si                  #14       Parawollastonite (CaSiO3) (part 3) 
AB20C2D14E2_mP78_14_a_10e_e_7e_e-001.{ABCDE}                         mP78      P2_{1}/c            H4_{4}        CuH20N2O14S2            #14       Tutton salt [Cu(NH4)2(SO4)2.6H2O] (part 3) 
A8B2CD12E2_mP100_14_8e_2e_ad_12e_2e-001.{ABCDE}                      mP100     P2_{1}/c            -             H8K2MnO12S2             #14       Manganese-leonite 110K [K2Mn(SO4)2.4H2O] (part 3) 
A_mC4_15_e-001.{A}                                                   mC4       C2/c                -             Ga                      #15       beta-Ga (obsolete) (part 3) 
ABC2_mC16_15_e_e_f-004.{ABC}                                         mC16      C2/c                F5_{a}        FeKS2                   #15       KFeS2 (part 3) 
A2B2C_mC20_15_ad_f_e-001.{ABC}                                       mC20      C2/c                -             Ag2O2Pb                 #15       Ag2PbO2 (part 3) 
A2B_mC24_15_2f_ce-001.{AB}                                           mC24      C2/c                -             O2Sb                    #15       Clinocervantite (beta-Sb2O4) (part 3) 
A2B5_mC28_15_f_e2f-003.{AB}                                          mC28      C2/c                -             Nb2O5                   #15       zeta-Nb2O5 (B-Nb2O5) (part 3) 
AB5C2_mC32_15_e_e2f_f-001.{ABC}                                      mC32      C2/c                -             NiSe5Ta2                #15       Ta2NiSe5 (part 3) 
AB5CD_mC32_15_e_e2f_e_b-001.{ABCD}                                   mC32      C2/c                S0_{6}        CaO5SiTi                #15       Titanite (CaTiSiO5) (part 3) 
A2BC4D2_mC36_15_f_e_2f_f-001.{ABCD}                                  mC36      C2/c                -             C2(H2O)O4Rb2            #15       Rb2C2O4.H2O (part 3) 
A7B2C2_mC44_15_e3f_f_f-001.{ABC}                                     mC44      C2/c                -             O7V2Zn2                 #15       alpha-Zn2V2O7 (part 3) 
AB4C6D_mC48_15_e_2f_3f_e-001.{ABCD}                                  mC48      C2/c                H4_{6}        CaH4O6S                 #15       Gypsum (CaSO4.2H2O) (part 3) 
AB4C4D4E_mC56_15_e_2f_2f_2f_a-001.{ABCDE}                            mC56      C2/c                H4_{22}       BaC4(H2O)4N4Ni          #15       BaNi(CN)4.4H2O (part 3) 
A5BC2_mC64_15_5f_f_2f-001.{ABC}                                      mC64      C2/c                -             O5SiY2                  #15       Y2SiO5 (RE2SiO5 X2) (part 3) 
AB5CD2_mC72_15_f_5f_f_2f-001.{ABCD}                                  mC72      C2/c                S5_{6}        AlO5(OH)Si2             #15       Pyrophyllite [AlSi2O5(OH)] (part 3) 
A2BC10D2E4_mC76_15_f_e_5f_f_2f-001.{ABCDE}                           mC76      C2/c                S5_{1}        Al3KO10(OH)2Si3         #15       Muscovite (KH2Al3Si3O12) (part 3) 
A2BCD12E3_mC76_15_f_e_b_6f_ef-001.{ABCDE}                            mC76      C2/c                -             Fe2MnNaO12P3            #15       Alluaudite [NaMnFe2(PO4)3] (part 3) 
A3B5C4D2_mC112_15_a3ef_5f_4f_2f-001.{ABCD}                           mC112     C2/c                S5_{5}        Mg3O5(OH)4Si2           #15       Chlorite [Mg3(Mg2Al)(Si3Al)O10(OH)8] (part 3) 
A2B4C2D17E6_mC124_15_f_2f_f_e8f_3f-001.{ABCDE}                       mC124     C2/c                -             BeHNaO8Si3              #15       Eudidymite (BeHNaO8Si3) (part 3) 
A2B3C9D3E_mC144_15_2f_bcdef_9f_3f_ae-001.{ABCDE}                     mC144     C2/c                -             (H2O)2Na2O9Si3Zr        #15       Catapleiite (Na2ZrSi3O9.2H2O) (part 3) 
A3B16C20D3_mC168_15_ef_8f_10f_ef-001.{ABCD}                          mC168     C2/c                H4_{20}       Cd3H16O20S3             #15       (CdSO4)3.8H2O (part 3) 
A8B2CD12E2_mC200_15_8f_2f_ce_2e11f_2f-001.{ABCDE}                    mC200     C2/c                -             H8K2MnO12S2             #15       Manganese-leonite 185K [K2Mn(SO4)2.4H2O] (part 3) 
ABC3_oP40_17_abcd_2e_abcd4e-001.{ABC}                                oP40      P222_{1}            -             NaNbO3                  #17       NaNbO3 (part 3) 
A2B_oP12_18_2c_c-001.{AB}                                            oP12      P2_{1}2_{1}2        -             O2Te                    #18       gamma-TeO2 (part 3) 
AB12C5D2_oP40_18_a_6c_b2c_c-001.{ABCD}                               oP40      P2_{1}2_{1}2        -             ClH12N5Zn2              #18       Diamminetriamidodizinc Chloride ([Zn2(NH3)2(NH2)3]Cl) (part 3) 
AB_oP16_19_2a_2a-001.{AB}                                            oP16      P2_{1}2_{1}2_{1}    -             NaP                     #19       NaP (part 3) 
A2B2C_oP20_19_2a_2a_a-001.{ABC}                                      oP20      P2_{1}2_{1}2_{1}    C31           H2O2Zn                  #19       W"ulfingite (epsilon-Zn(OH)2) (part 3) 
AB4C_oP24_19_a_4a_a-001.{ABC}                                        oP24      P2_{1}2_{1}2_{1}    -             AlCl4Na                 #19       NaAlCl4 (part 3) 
A6BC4D_oP48_19_6a_a_4a_a-001.{ABCD}                                  oP48      P2_{1}2_{1}2_{1}    -             H6NO4P                  #19       Ferroelectric NH4H2PO4 (part 3) 
AB2C_oP80_19_5a_10a_5a-001.{ABC}                                     oP80      P2_{1}2_{1}2_{1}    -             CH2O                    #19       beta-Arabinose [(CH2O)20] (part 3) 
A14BC11D_oP108_19_14a_a_11a_a-001.{ABCD}                             oP108     P2_{1}2_{1}2_{1}    H4_{12}       H14NiO11S               #19       Morenosite (NiSO4.7H2O) (part 3) 
AB3_oC16_20_a_bc-001.{AB}                                            oC16      C222_{1}            D0_{7}        CrO3                    #20       D07 (CrO3) (obsolete) (part 3) 
AB4C_oC24_20_b_2c_a-001.ABC/CAB/BCA                                  oC24      C222_{1}            -             AlO4P                   #20       AlPO4 ``low cristobalite type'' (part 3) 
AB5C2_oC32_20_b_a2bc_c-001.{ABC}                                     oC32      C222_{1}            K3_{3}        AlF5Tl2                 #20       Tl2AlF5 (part 3) 
AB2_oC6_21_a_k-002.{AB}                                              oC6       C222                -             HoSb2                   #21       HoSb2 (part 3) 
A2B7C2_oF88_22_k_bdefghij_k-001.{ABC}                                oF88      F222                -             Cd2O7Re2                #22       Predicted Phase IV Cd2Re2O7 (part 3) 
AB6_oP28_29_a_6a-001.{AB}                                            oP28      Pca2_{1}            -             HgN6                    #29       Mercury (II) Azide [Hg(N3)2] (part 3) 
ABC30DE20F2_oP220_29_a_a_30a_a_20a_2a-001.{ABCDEF}                   oP220     Pca2_{1}            -             AlCH30NO20S2            #29       Low-Temperature (NH3CH3)Al(SO4)2.12H2O (part 3) 
A4B7C_oP24_31_2b_a3b_a-001.{ABC}                                     oP24      Pmn2_{1}            -             B4O7Sr                  #31       B4SrO7 (part 3) 
A2B6CD8_oP34_31_2a_2a2b_a_4a2b-001.{ABCD}                            oP34      Pmn2_{1}            H4_{11}       Cl2(H2O)6MgO8           #31       Mg(ClO4)2.6H2O (part 3) 
A13B4_oP102_31_17a11b_8a2b-001.{AB}                                  oP102     Pmn2_{1}            -             Al13Co4                 #31       Orthorhombic Co4Al13 (part 3) equivalent to 491,492 [HTQC]
A17B47_oP128_32_a8c_a23c-001.{AB}                                    oP128     Pba2                -             Mo17O47                 #32       Mo17O47 (part 3) 
ABC3_oP20_33_a_a_3a-001.{ABC}                                        oP20      Pna2_{1}            -             ILiO3                   #33       gamma-LiIO3 (part 3) 
A2B7C2_oP44_33_2a_7a_2a-001.{ABC}                                    oP44      Pna2_{1}            -             Gd2O7Si2                #33       Possible delta-Gd2Si2O7 (part 3) 
A4BCD6_oP48_33_4a_a_a_6a-001.{ABCD}                                  oP48      Pna2_{1}            -             B4CsFO6                 #33       CsB4O6F (part 3) 
A2BC4_oP84_33_6a_3a_12a-001.{ABC}                                    oP84      Pna2_{1}            -             B2CaO4                  #33       CaB2O4 (III) (part 3) 
A2BC2_oP10_34_c_a_c-001.{ABC}                                        oP10      Pnn2                -             F2-xHxMnOx              #34       MnF2-x(OH)x (part 3) 
AB2_oC12_36_a_2a-002.{AB}                                            oC12      Cmc2_{1}            -             MoP2                    #36       MoP2 (part 3) 
A2BC2_oC20_36_b_a_b-001.{ABC}                                        oC20      Cmc2_{1}            -             N2OSi2                  #36       Si2N2O (part 3) 
A2BC5_oC32_36_b_a_a2b-001.{ABC}                                      oC32      Cmc2_{1}            -             Bi2GeO5                 #36       Bi2GeO5 (part 3) 
A4B7C2D2_oC60_36_2b_a3b_2a_b-001.{ABCD}                              oC60      Cmc2_{1}            S4_{6}        Be4O7(OH)2Si2           #36       Bertrandite (Be4Si2O7(OH)2) (part 3) 
A3B2_oC80_36_4a4b_2a3b-001.{AB}                                      oC80      Cmc2_{1}            -             Ni3Si2                  #36       Ni3Si2 (part 3) 
ABC3_oC80_36_2ab_2ab_2a5b-001.{ABC}                                  oC80      Cmc2_{1}            -             KNO3                    #36       alpha-Potassium Nitrate (KNO3) II (part 3) 
ABC6D15_oC46_38_b_b_2a2d_2ab4d2e-001.{ABCD}                          oC46      Amm2                -             FNaNb6O15               #38       NaNb6O15F (part 3) 
AB3_oC16_40_b_a2b-002.{AB}                                           oC16      Ama2                -             CrO3                    #40       Orthorhombic CrO3 (part 3) 
A2B7C2_oC88_40_abc_2b6c_a3b-001.{ABC}                                oC88      Ama2                -             Mo2O7Rb2                #40       Rb2Mo2O7 (part 3) 
A5B8CD12_oC104_41_a2b_4b_a_6b-001.{ABCD}                             oC104     Aea2                K3_{5}        B5H8KO12                #41       Santite (KB5O8.4H2O) (part 3) 
A2B_oF24_43_b_a-001.{AB}                                             oF24      Fdd2                -             Cs2Se                   #43       Cs2Se (part 3) 
A3B2_oF40_43_ab_b-001.{AB}                                           oF40      Fdd2                -             Al3Zr2                  #43       Zr2Al3 (part 3) equivalent to 365,366 [HTQC]
A2BC4D_oF64_43_b_a_2b_a-001.{ABCD}                                   oF64      Fdd2                -             H2KO4P                  #43       Archerite (KH2PO4) (part 3) 
A2B7C2_oF88_43_b_a3b_b-001.{ABC}                                     oF88      Fdd2                -             Cu2O7V2                 #43       Blossite (alpha-Cu2V2O7) (part 3) 
A2B4C2D12E3_oF184_43_b_2b_b_6b_ab-001.{ABCDE}                        oF184     Fdd2                S6_{10}       Al2H4Na2O12Si3          #43       Natrolite (Na2Al2Si3O10.2H2O) (part 3) 
ABC2_oI8_44_a_a_c-001.{ABC}                                          oI8       Imm2                F5_{5}        NaNO2                   #44       Ferroelectric NaNO2 (part 3) 
ABC2_oI8_44_a_a_d-001.{ABC}                                          oI8       Imm2                F5_{12}       AgNO2                   #44       AgNO2 (part 3) 
A2B5CD2_oI40_44_2c_abcde_d_e-001.{ABCD}                              oI40      Imm2                S2_{2}        H2O5SiZn2               #44       Hemimorphite (Zn4Si2O7(OH)2.H2O) (part 3) 
AB_oI48_44_6d_ab2cde-001.{AB}                                        oI48      Imm2                B30           MgZn                    #44       B30 (MgZn?) (part 3) 
A2B17C6_oI100_46_ab_b8c_3c-001.{ABC}                                 oI100     Ima2                -             Nb2O17Zr6               #46       Nb2Zr6O17 (part 3) 
AB2C_oP8_51_e_be_f-001.{ABC}                                         oP8       Pmma                -             Bi2Ni3S2                #51       Parkerite (Ni3Bi2S2) (part 3) 
ABC6D15_oP46_51_f_d_2e2i_aef4i2j-001.{ABCD}                          oP46      Pmma                -             FLiNb6O15               #51       LiNb6O15F (part 3) 
A3B12CDE6_oP276_52_d4e_18e_ce_de_2d8e-001.{ABCDE}                    oP276     Pnna                -             Cl3H12KMgO6             #52       Carnallite [Mg(H2O)6KCl3] (part 3) 
A2BC_oP16_53_eh_ab_g-001.{ABC}                                       oP16      Pmna                F5_{8}        F2H5N                   #53       NH4HF2 (part 3) 
A2BC4D2_oP18_53_h_a_i_e-001.{ABCD}                                   oP18      Pmna                C45           Cl2CuH4O2               #53       Eriochalcite (CuCl2.2H2O) (part 3) 
A2B3C8_oP26_55_h_ag_2g2h-001.{ABC}                                   oP26      Pbam                -             Nb2Pd3Se8               #55       Nb2Pd3Se8 (part 3) 
A4BCD2_oP32_55_ghi_f_e_gh-001.{ABCD}                                 oP32      Pbam                E3_{4}        Cl4(H2O)HgK2            #55       K2HgCl4.H2O (part 3) 
AB2C5_oP32_55_g_fh_eghi-001.{ABC}                                    oP32      Pbam                -             HoMn2O5                 #55       HoMn2O5 (part 3) 
A8B11_oP38_55_g3h_a3g2h-001.{AB}                                     oP38      Pbam                -             B8Ru11                  #55       Ru11B8 (part 3) equivalent to 1086.AB,1086.BA [HTQC]
A10B3C4_oP68_55_2e2fgh2i_adef_2e2f-001.{ABC}                         oP68      Pbam                -             O10Ru3Sr4               #55       Orthorhombic Sr4Ru3O10 (part 3)
A2BC4_oP56_56_2e_e_4e-001.{ABC}                                      oP56      Pccn                -             B2CaO4                  #56       Calciborite (CaB2O4 II) (part 3) 
A3B_oP16_57_a2d_d-001.{AB}                                           oP16      Pbcm                D0_{10}       O3W                     #57       D010 (WO3) (em obsolete) (part 3) 
A4BC_oP24_57_cde_d_a-001.{ABC}                                       oP24      Pbcm                -             O4SrU                   #57       SrUO4 (part 3) 
ABC3_oP40_57_cd_e_cd2e-001.{ABC}                                     oP40      Pbcm                -             NaNbO3                  #57       Lueshite (NaNbO3) (part 3) 
AB_oP8_58_g_g-001.{AB}                                               oP8       Pnnm                -             InS                     #58       InS (part 3) 
A2B3C6_oP22_58_g_af_gh-001.{ABC}                                     oP22      Pnnm                -             B2Mg3O6                 #58       Kotoite (Mg3(BO3)2) (part 3) 
A4B3_oP28_58_4g_3g-001.{AB}                                          oP28      Pnnm                -             In4Se3                  #58       In4Se3 (part 3) 
ABC5D2_oP36_58_g_g_3gh_eg-001.{ABCD}                                 oP36      Pnnm                H2_{7}        AsHO5Zn2                #58       Adamite [Zn2(AsO4)(OH)] (part 3) 
A2B7C24D8_oP82_58_g_ae2f_2g5h_2h-001.{ABCD}                          oP82      Pnnm                -             H2Mg7O24Si8             #58       Protoanthophyllite (H2Mg7Si8O24) (part 3) 
A4B2C3_oP18_59_ef_ab_af-001.{ABC}                                    oP18      Pmmn                G0_{11}       H4N2O3                  #59       NH4NO3 IV (part 3) 
A2BC4_oP28_60_d_c_2d-001.{ABC}                                       oP28      Pbcn                E3_{2}        B2CaO4                  #60       CaB2O4 I (part 3) 
AB2C6_oP36_60_c_d_3d-001.{ABC}                                       oP36      Pbcn                E5_{1}        FeNb2O6                 #60       Columbite (FeNb2O6) (part 3) 
A5B12_oP68_60_c2d_6d-001.{AB}                                        oP68      Pbcn                -             Cr5O12                  #60       Cr5O12 (part 3) 
A2B_oP24_61_2c_c-002.{AB}                                            oP24      Pbca                C52           O2Te                    #61       Tellurite (beta-TeO2) (part 3) 
ABC_oP24_61_c_c_c-001.{ABC}                                          oP24      Pbca                -             CClO                    #61       COCl (part 3) 
A2B4C_oP28_61_c_2c_a-001.{ABC}                                       oP28      Pbca                -             Ca2O4Ru                 #61       Ca2RuO4 (part 3) 
AB2CD4_oP64_61_c_2c_c_4c-001.{ABCD}                                  oP64      Pbca                G7_{2}        BBe2HO4                 #61       Hambergite [Be2BO3(OH)] (part 3) 
A7BCD_oP80_61_7c_c_c_c-001.{ABCD}                                    oP80      Pbca                -             Cl7OPTi                 #61       (TiCl4.POCl3)2 (part 3) 
AB3C_oP80_61_2c_6c_2c-001.{ABC}                                      oP80      Pbca                S4_{3}        MgO3Si                  #61       Enstatite (MgSiO3) (part 3) 
A2B_oP12_62_2c_c-010.{AB}                                            oP12      Pnma                C53           Br2Sr                   #62       C53 (SrBr2) (obsolete) (part 3) 
A3B_oP16_62_3c_c-002.{AB}                                            oP16      Pnma                D0_{16}       I3(NH4)                 #62       NH4I3 (part 3) 
ABC2_oP16_62_c_c_2c-002.{ABC}                                        oP16      Pnma                E0_{2}        AlHO2                   #62       Diaspore (AlOOH) (part 3) 
ABCD_oP16_62_c_c_c_c-001.{ABCD}                                      oP16      Pnma                F5_{14}       ClBrI(NH4)              #62       NH4ClBrI (part 3) 
AB3C_oP20_62_c_3c_c-001.{ABC}                                        oP20      Pnma                E2_{4}        CdCl3(NH4)              #62       NH4CdCl3 (part 3) 
ABC3_oP20_62_c_c_cd-001.{ABC}                                        oP20      Pnma                -             KNO3                    #62       alpha-Potassium Nitrate (KNO3) I (part 3) 
ABC3_oP20_62_c_c_cd-002.{ABC}                                        oP20      Pnma                G0_{10}       N(NH4)O3                #62       NH4NO3 III (part 3) 
ABC3_oP20_62_c_c_cd-003.{ABC}                                        oP20      Pnma                G0_{2}        CaCO3                   #62       Aragonite (CaCO3) (part 3) 
A2B_oP24_62_4c_2c-001.{AB}                                           oP24      Pnma                -             Cs2Sb                   #62       Cs2Sb (part 3) 
AB4C_oP24_62_a_2cd_c-001.{ABC}                                       oP24      Pnma                -             CuO4S                   #62       Chalcocyanite (CuSO4) (part 3) 
A2B4C_oP28_62_2c_2cd_c-001.{ABC}                                     oP28      Pnma                H1_{6}        K2O4S                   #62       Arcanite (K2SO4) (part 3) 
A5BC_oP28_62_3cd_c_c-001.{ABC}                                       oP28      Pnma                -             O5SV                    #62       VOSO4 (part 3) 
AB4C2_oP28_62_c_4c_2c-001.{ABC}                                      oP28      Pnma                E3_{3}        FeS4Sb2                 #62       Berthierite (FeSb2S4) (part 3) 
AB6_oP28_62_c_6c-001.{AB}                                            oP28      Pnma                -             CuN6                    #62       Copper (II) Azide [Cu(N3)2] (part 3) 
A2B5C_oP32_62_bc_3cd_c-001.{ABC}                                     oP32      Pnma                S0_{3}        Al2O5Si                 #62       Sillimanite (Al2SiO5) (part 3) 
A3B_oP32_62_ab4c_2c-001.{AB}                                         oP32      Pnma                -             O3W                     #62       Original beta-WO3 (em obsolete) (part 3) 
A4BC2D_oP32_62_2cd_b_2c_a-001.{ABCD}                                 oP32      Pnma                E3_{5}        Cl4(H2O)K2Sn            #62       K2SnCl4.H2O (part 3) 
A4BC2D_oP32_62_2cd_c_d_c-001.{ABCD}                                  oP32      Pnma                -             Cl4(H2O)K2Sn            #62       K2SnCl4.H2O (part 3) 
A2B2C4D_oP36_62_d_d_2cd_c-001.{ABCD}                                 oP36      Pnma                S0_{5}        Al2F2O4Si               #62       Topaz (Al2SiO4F2, S05) (part 3) 
AB2C3D3_oP36_62_c_ac_cd_cd-001.{ABCD}                                oP36      Pnma                -             ClCu2H3O3               #62       Atacamite (Cu2(OH)3Cl) (part 3) 
AB6C2_oP36_62_c_2c2d_d-001.{ABC}                                     oP36      Pnma                -             CaO6Ta2                 #62       Rynersonite (Orthorhombic CaTa2O6) (part 3) 
A2B3C4D_oP40_62_d_cd_2cd_c-001.{ABCD}                                oP40      Pnma                S0_{7}        F2Mg3O4Si               #62       Norbergite [Mg(F,OH)2.Mg2SiO4] (part 3) 
A2B6C3_oP44_62_2c_2c2d_3c-001.{ABC}                                  oP44      Pnma                K5_{1}        K2O6S3                  #62       K2S3O6 (part 3) 
A7B2C2_oP44_62_3c2d_2c_d-001.{ABC}                                   oP44      Pnma                -             O7Si2Y2                 #62       Possible delta-Y2Si2O7 (part 3) 
A8BCD_oP44_62_4c2d_c_c_c-001.{ABCD}                                  oP44      Pnma                -             Cl8OPSb                 #62       SbCl5.POCl3 (part 3) 
A2BC8D2_oP52_62_d_c_2c3d_d-001.{ABCD}                                oP52      Pnma                S6_{3}        B2CaO8Si2               #62       Danburite (CaB2Si2O8) (part 3) 
A4B3_oP56_62_8c_6c-001.{AB}                                          oP56      Pnma                -             Mo4P3                   #62       Mo4P3 (part 3) 
A3B15C5D_oP96_62_cd_3c6d_3cd_c-001.{ABCD}                            oP96      Pnma                J1_{8}        Cl3H15N5Rh              #62       RhCl2(NH3)5Cl (part 3) 
A8B4C4DE8F2_oP108_62_4c2d_2d_2cd_c_4c2d_d-001.{ABCDEF}               oP108     Pnma                F2_{1}        C8H4K4MoN8O2            #62       K4[Mo(CN)8].2H2O (part 3) 
A4B3_oP112_62_8c4d_4c4d-001.{AB}                                     oP112     Pnma                -             P4Se3                   #62       P4Se3 (part 3) 
ABCD8E3_oP112_62_d_2c_d_4c6d_3d-001.{ABCDE}                          oP112     Pnma                S4_{7}        BeHNaO8Si3              #62       Epididymite (BeHNaO8Si3) (part 3) 
A2B5C22D2E8_oP156_62_d_c2d_2c10d_2c_4d-001.{ABCDE}                   oP156     Pnma                S4_{4}        Fe2Mg5O22(OH)2Si8       #62       Anthophyllite (Mg5Fe2Si8O22(OH)2) (part 3) 
AB22C23D2E2_oP200_62_c_11d_3c10d_d_d-001.{ABCDE}                     oP200     Pnma                -             CaH22O23P2U2            #62       Autunite Ca[(UO2)(PO4)]2(H2O)11 (part 3) 
AB3_oC16_63_c_3c-001.{AB}                                            oC16      Cmcm                D0_{d}        AsMn3                   #63       Mn3As (part 3) 
AB2C2_oC20_63_c_f_2c-001.{ABC}                                       oC20      Cmcm                E0_{4}        FeHO2                   #63       Lepidocrocite (gamma-FeO(OH)) (part 3) 
AB2CD_oC20_63_b_f_c_c-001.{ABCD}                                     oC20      Cmcm                -             CFe2SiTh                #63       ThFe2SiC (part 3) 
A5B_oC24_63_c2f_c-001.{AB}                                           oC24      Cmcm                -             Te5Zr                   #63       ZrTe5 (part 3) 
A_oC24_63_3f-001.{A}                                                 oC24      Cmcm                -             Si24                    #63       Si24 Clathrate (part 3) 
A5B3_oC32_63_cfg_ce-001.{AB}                                         oC32      Cmcm                -             Pd5Pu3                  #63       Pd5Pu3 (part 3) equivalent to 661.AB,661.BA [HTQC]
AB5C2_oC32_63_c_c2f_f-004.{ABC}                                      oC32      Cmcm                -             NiS5Ta2                 #63       Ta2NiS5 (part 3) 
A2B2C6DE2_oC52_63_g_e_fh_c_f-001.{ABCDE}                             oC52      Cmcm                -             Br2Cu2O6PbSe2           #63       Cu2Pb(SeO3)2Br2 (part 3) 
A4BC12D2_oC76_63_eg_c_f3gh_g-001.{ABCD}                              oC76      Cmcm                S0_{4}        Al4FeO12Si2             #63       S04 (Staurolite, Fe(OH)2Al4Si2O10) (obsolete) (part 3) 
A10B3C4_oC68_64_2dfg_ad_2d-001.{ABC}                                 oC68      Cmce                -             O10Ru3Sr4               #64       Base-centered orthorhombic Sr4Ru3O10 (part 3) 
A2B2C7_oC88_64_ef_df_3f2g-001.{ABC}                                  oC88      Cmce                -             Mo2Na2O7                #64       Na2Mo2O7 (part 3) 
A3B8_oC22_65_ag_bd2gh-001.{AB}                                       oC22      Cmmm                -             FNb3O7                  #65       Nb3O7F (part 3) 
A2B8CD2_oC26_65_h_r_a_i-001.{ABCD}                                   oC26      Cmmm                E1_{3}        Cl2H6MgN2               #65       Mg(NH3)2Cl2 (part 3) 
A2BC2D_oC24_67_m_a_n_g-001.{ABCD}                                    oC24      Cmme                F5_{7}        H2(NH4)O2P              #67       NH4H2PO2 (part 3) 
A2B4C_oF56_70_g_h_a-001.{ABC}                                        oF56      Fddd                H1_{7}        Na2O4S                  #70       Thenardite [Na2SO4 (V)] (part 3) 
ABC2_oI16_71_g_i_eh-001.{ABC}                                        oI16      Immm                -             CsFeS2                  #71       CsFeS2 (100K) (part 3) 
AB6C3_oI20_71_a_in_cj-001.{ABC}                                      oI20      Immm                -             AlF6Na3                 #71       High-Temperature Cryolite (Na3AlF6) (part 3) 
A2B5_oI28_72_j_bfj-001.{AB}                                          oI28      Ibam                D8_{g}        Ga2Mg5                  #72       Ga2Mg5 (part 3) 
ABC4D_oI28_74_a_d_hi_e-001.{ABCD}                                    oI28      Imma                -             CuLiO4V                 #74       LiCuVO4 (part 3) 
A2B6C2D_oI44_74_h_ij_i_e-001.{ABCD}                                  oI44      Imma                E1_{2}        Cl2H6N2Zn               #74       Zn(NH3)2Cl2 (part 3) 
A4B2C3_tP72_77_8d_ab2c2d_6d-001.{ABC}                                tP72      P4_{2}              -             H4N2O3                  #77       Gwihabaite [NH4NO3 (V)] (part 3) 
AB5C_tP14_85_c_cg_b-001.{ABC}                                        tP14      P4/n                -             MoO5P                   #85       MoPO5 (part 3) 
A3B6CD_tP44_85_bcg_3g_ac_e-001.{ABCD}                                tP44      P4/n                E2_{6}        (Br,Cl)3(H2O)6KMg       #85       Bromocarnallite (KMg(H2O)6(Cl,Br)3) (part 3) 
A2BC_tP32_86_2g_g_g-001.{ABC}                                        tP32      P4_{2}/n            E1_{4}        Cl2NP                   #86       PNCl2 (part 3) 
AB6C_tP32_86_d_3g_c-001.{ABC}                                        tP32      P4_{2}/n            J1_{11}       Na(OH)6Sb               #86       NaSb(OH)6 (part 3) 
ABC3_tP40_86_g_g_3g-001.{ABC}                                        tP40      P4_{2}/n            -             ILiO3                   #86       beta-LiIO3 (part 3) 
A4B11C2_tP68_86_2g_ab5g_g-001.{ABC}                                  tP68      P4_{2}/n            -             Nd4O11Re2               #86       Nd4Re2O11 (part 3) 
AB6C2D_tI20_87_a_eh_d_b-001.ABCD/DACB/CABD/                          tI20      I4/m                -             NiO6Sr2W                #87       Sr2NiWO6 (part 3) 
DCBA/BCAD/DBAC/CDAB/ADBC/BADC/CBDA/ACDB/BDCA
AB4C24D12_tI82_87_a_h_2h2i_hi-001.{ABCD}                             tI82      I4/m                S6_{4}        ClNa4O24(Al3Si9)        #87       Marialite Scapolite [Na4Cl(AlSi3)3O24] (part 3) 
A9B4C20_tI132_88_a2f_f_5f-001.{ABC}                                  tI132     I4_{1}/a            -             Ge9Na4O20               #88       Na4Ge9O20 (part 3) 
AB2C3D2_tP16_90_c_f_ce_e-001.{ABCD}                                  tP16      P42_{1}2            G7_{5}        CCl2O3Pb2               #90       G75 (PbCO3.PbCl2, Phosgenite) (obsolete) (part 3) 
A12BC10D_tP96_92_6b_a_5b_a-001.{ABCD}                                tP96      P4_{1}2_{1}2        H4_{5}        H12NiO10S               #92       Retgersite (alpha-NiSO4.6H2O) (part 3) 
A2B7C2_tI44_98_f_bcde_f-001.{ABC}                                    tI44      I4_{1}22            -             Cd2O7Re2                #98       Phase III Cd2Re2O7 (part 3) 
ABC2_tP8_100_b_a_c-001.{ABC}                                         tP8       P4bm                F5_{4}        Cl(NH4)O2               #100      F54 (NH4ClO2) (obsolete) (part 3) 
ABC3_tP10_100_b_a_bc-001.{ABC}                                       tP10      P4bm                G0_{9}        N(NH4)O3                #100      NH4NO3 II (part 3) 
A6B2C_tP72_103_abc5d_2d_abc-001.{ABC}                                tP72      P4cc                -             O6Se2V                  #103      VSe2O6 (part 3) 
A2BC4_tP7_111_f_a_n-001.{ABC}                                        tP7       P-42m               E3_{1}        Ag2HgI4                 #111      E31 (beta-Ag2HgI4) (obsolete) (part 3) 
AB4CD2_tP16_113_c_f_a_e-001.{ABCD}                                   tP16      P-42_{1}m           -             ClH4NO2                 #113      Ammonium Chlorite (NH4ClO2) (part 3) 
A2B12C4D4E_tP46_114_d_3e_e_e_a-001.{ABCDE}                           tP46      P-42_{1}c           H4_{17}       Ag2H12N4O4S             #114      Ag2SO4.4NH3 (part 3) 
A19B15_tP68_114_bc4e_ac3e-001.{AB}                                   tP68      P-42_{1}c           -             C19Sc15                 #114      C19Sc15 (part 3) equivalent to 1217.AB,1217.BA,1346.AB,1346.BA [HTQC]
ABC2_tP4_115_a_c_g-001.{ABC}                                         tP4       P-4m2               F6_{1}        CuFeS2                  #115      F61 (Chalcopyrite, CuFeS2) (obsolete) (part 3) 
AB2C_tI8_119_c_e_a-001.{ABC}                                         tI8       I-4m2               -             FeS2Tl                  #119      Tetragonal TlFeS2 (part 3) 
A2B7C2_tI44_119_i_bdefgh_i-001.{ABC}                                 tI44      I-4m2               -             Cd2O7Re2                #119      Phase II Cd2Re2O7 (part 3)
AB8C8D_tI72_120_c_2i_2i_b-001.{ABCD}                                 tI72      I-4c2               H4_{3}        BeH8O8S                 #120      BeSO4.4H2O (part 3) 
AB3_tI32_121_g_f2i-001.{AB}                                          tI32      I-42m               -             SV3                     #121      alpha-V3S (part 3) equivalent to 649.AB,649.BA [HTQC]
A2B2C6D_tI44_121_i_i_ij_c-001.{ABCD}                                 tI44      I-42m               -             B2Cu2O6Sr               #121      SrCu2(BO3)2 (part 3) 
A2BC2_tI40_122_e_d_e-001.{ABC}                                       tI40      I-42d               F1_{1}        C2HgN2                  #122      Mercury Cyanide [Hg(CN)2] (part 3) 
A4BC4D_tI40_122_e_b_e_a-001.{ABCD}                                   tI40      I-42d               H2_{2}        H2KO4P                  #122      KH2PO4 (part 3) 
AB2_tI48_122_cd_2e-001.{AB}                                          tI48      I-42d               -             NaS2                    #122      NaS2 (part 3) 
A8BC4D_tI56_122_2e_b_e_a-001.{ABCD}                                  tI56      I-42d               -             H6NO4P                  #122      NH4H2PO4 (part 3) 
AB_tP2_123_a_d-002.AB                                                tP2       P4/mmm              L2_{a}        CuTi                    #123      delta-CuTi (part 3) 
A3BC_tP5_123_cg_a_d-001.{ABC}                                        tP5       P4/mmm              E2_{5}        Cl3Hg(NH4)              #123      NH4HgCl3 (part 3) 
AB4C_tP6_123_d_eh_a-001.{ABC}                                        tP6       P4/mmm              H0_{8}        AlF4Tl                  #123      TlAlF4 (part 3) 
A4B2C_tP7_123_j_e_a-001.{ABC}                                        tP7       P4/mmm              H1_{5}        Cl4K2Pt                 #123      K2PtCl4 (part 3) 
A8B2C_tP11_123_r_f_a-001.{ABC}                                       tP11      P4/mmm              E6_{1}        Sr(OH)2(H2O)8           #123      E61 (Sr(OH)2(H2O)8) (obsolete) (part 3) 
A8B2C_tP11_123_r_h_a-001.{ABC}                                       tP11      P4/mmm              E6_{2}        (H2O)8SrO2              #123      E62 [SrO2(H2O)8] (empossibly obsolete) (part 3) 
AB8C2_tP22_124_a_n_h-001.{ABC}                                       tP22      P4/mcc              -             Ca(H2O)8O2              #124      CaO2(H2O)8 (part 3) 
ABC4D2E8_tP32_126_a_b_h_e_k-001.{ABCDE}                              tP32      P4/nnc              J1_{9}        AgCoN4(NH3)2O8          #126      Ag[Co(NH3)2(NO2)4] (part 3) 
A4B10C2D34E4F9_tP252_126_k_ce2k_f_h8k_k_d2k-001.{ABCDEF}             tP252     P4/nnc              S2_{3}        Ag4Ca10Mg2O34(OH)4Si9   #126      Vesuvianite (Ca10Al4(Mg,Fe)2Si9O34(OH)4) (part 3) 
A2BC4D_tP16_127_h_d_i_a-001.{ABCD}                                   tP16      P4/mbm              H4_{9}        Cl2(H2O)(NH3)4Pd        #127      Pd(NH3)4Cl2.H2O (part 3) 
AB2C3D2_tP32_127_g_eh_gk_k-001.{ABCD}                                tP32      P4/mbm              -             CCl2O3Pb2               #127      Phosgenite [Pb2Cl2(CO3)] (part 3) 
A3B14C5_tP44_128_ac_ehi_bg-001.{ABC}                                 tP44      P4/mnc              K7_{5}        Al3F14Na5               #128      Chiolite (Na5Al3F14) (part 3) 
A4BC16DE28F8_tP116_128_h_a_2i_b_g3i_i-001.{ABCDEF}                   tP116     P4/mnc              S5_{2}        Ca4FH16KO28Si8          #128      Apophyllite (KCa4Si8O20F.8H2O) (part 3) 
ABCD_tP8_129_b_c_a_c-001.{ABCD}                                      tP8       P4/nmm              -             AgLaOS                  #129      LaOAgS (part 3) 
AB4C_tP12_129_c_i_a-001.{ABC}                                        tP12      P4/nmm              B25           BrH4N                   #129      NH4Br (part 3) 
AB4C6DE_tP26_129_c_j_2ci_a_c-001.{ABCDE}                             tP26      P4/nmm              H5_{10}       Ca(H2O)6O12P2U2         #129      Meta-autunite (I) [Ca(UO2)2(PO4)2.6H2O] (part 3) 
A18B10C_tP116_130_2c4g_2c2g_a-001.{ABC}                              tP116     P4/ncc              -             H18O10Sr                #130      Sr(OH)2(H2O)8 (part 3) 
A4B2C_tP14_136_i_g_b-001.{ABC}                                       tP14      P4_{2}/mnm          -             Fe4Si2Zr                #136      ZrFe4Si2 (part 3) 
A2B3_tP20_136_j_dfg-001.{AB}                                         tP20      P4_{2}/mnm          -             Al2Zr3                  #136      Zr3Al2 (part 3) equivalent to 361,362 [HTQC]
A4BC4D2E2_tP26_136_fg_a_j_d_e-001.{ABCDE}                            tP26      P4_{2}/mnm          H4_{1}        Cl4CuH4K2O2             #136      K2CuCl4.2H2O (part 3) 
AB14C2_tP68_136_f_ce2j2k_fg-001.{ABC}                                tP68      P4_{2}/mnm          -             BFe14Nd2                #136      Nd2Fe14B (part 3) 
A2B3_tI10_139_e_ae-001.{AB}                                          tI10      I4/mmm              -             Au2Nb3                  #139      Au2Nb3 (part 3) 
AB2C2_tI10_139_a_d_e-003.{ABC}                                       tI10      I4/mmm              H5_{9}        CaP2U2                  #139      H59 [Autunite, Ca(UO2)2(PO4)2.10(1/2)H2O] (obsolete) (part 3) 
A4B2C2D_tI18_139_h_d_e_a-001.{ABCD}                                  tI18      I4/mmm              J1_{5}        Cl4K2O2Os               #139      K2OsO2Cl4 (part 3) 
AB3C_tI20_139_ab_eh_d-001.{ABC}                                      tI20      I4/mmm              K7_{6}        AuCl3Cs                 #139      AuCsCl3 (part 3) 
A10B4C3_tI34_139_c2eg_2e_ae-001.{ABC}                                tI34      I4/mmm              -             O10Sr4Ti3               #139      Sr4Ti3O10 (part 3) 
A6B2C3D_tI168_139_egikl2m_ejn_bh2n_acf-001.{ABCD}                    tI168     I4/mmm              J3_{1}        Cl6(H2O)2K3Tl           #139      K3TlCl6.2H2O (part 3) 
A2BC_tI16_140_h_d_a-001.{ABC}                                        tI16      I4/mcm              F5_{2}        F2HK                    #140      KHF2 (part 3) 
AB6_tI28_140_a_hk-001.{AB}                                           tI28      I4/mcm              D2_{c}        MnU6                    #140      U6Mn (part 3) 
A5BC2_tI32_140_bl_a_h-001.{ABC}                                      tI32      I4/mcm              K3_{4}        Br5(NH4)Pb2             #140      NH4Pb2Br5 (part 3) 
A5BC3_tI36_140_cl_b_ah-001.{ABC}                                     tI36      I4/mcm              K3_{1}        Cl5CoCs3                #140      Cs3CoCl5 (part 3) 
A31B20_tI204_140_b2gh3m_ac2fh3l-001.{AB}                             tI204     I4/mcm              -             Pu31Rh20                #140      Pu31Rh20 (part 3) 
AB11_tI48_141_a_bdi-001.{AB}                                         tI48      I4_{1}/amd          -             BaCd11                  #141      BaCd11 (part 3) 
A2B3_tI160_142_deg_3g-001.{AB}                                       tI160     I4_{1}/acd          -             As2Cd3                  #142      Cd3As2 (part 3) 
A2B2C3D12E4_tI184_142_f_f_be_3g_g-001.{ABCDE}                        tI184     I4_{1}/acd          S6_{1}        Al(H2O)NaO6Si2          #142      Analcime (NaAlSi2O6.H2O) (part 3) 
AB3C9D_hP28_143_2a_2d_6d_bc-001.{ABCD}                               hP28      P3                  -             BLa3O9W                 #143      La3BWO9 (P3) (part 3) 
AB3C_hP45_144_3a_9a_3a-001.{ABC}                                     hP45      P3_{1}              -             NO3Rb                   #144      RbNO3 (IV) (part 3) 
A2B3C_hP12_147_abd_g_d-001.{ABC}                                     hP12      P-3                 G3_{2}        Na2O3S                  #147      Na2SO3 (part 3) 
A2BCD6_hR10_148_c_a_b_f-001.{ABCD}                                   hR10      R-3                 G1_{1}        C2CaMgO6                #148      Dolomite [MgCa(CO3)2] (part 3) 
A6B6CD_hR14_148_f_f_b_a-001.ABCD/CABD/BCAD/                          hR14      R-3                 I6_{1}        Cl6(H2O)6NiSn           #148      Ni(H2O)6SnCl6 (part 3) 
DBAC/ADBC/CDAB/BACD/ACBD/CBAD/BDAC/DABC/DCAB
A6B2C6D_hR15_148_f_c_f_a-001.{ABCD}                                  hR15      R-3                 H6_{2}        H6K2O6Sn                #148      K2Sn(OH)6 (part 3) 
ABC8D2_hP12_150_b_a_dg_d-001.{ABCD}                                  hP12      P321                H3_{2}        AlKO8S2                 #150      Steklite [KAl(SO4)2] (part 3) 
A2B9C3_hP14_150_d_eg_ad-001.{ABC}                                    hP14      P321                K7_{3}        As2Cl9Cs3               #150      Cs3As2Cl9 (part 3) 
A2B12C6D_hP21_150_d_2g_ef_a-001.{ABCD}                               hP21      P321                -             Cl2H12O6Sr              #150      SrCl2.(H2O)6 (part 3) 
AB2CD6_hP30_150_e_c2d_f_3g-001.{ABCD}                                hP30      P321                -             BaC2CaO6                #150      Paralstonite (BaCa(CO3)2) (part 3) 
AB3C_hP30_150_ef_3g_c2d-001.{ABC}                                    hP30      P321                K1_{1}        KO3S                    #150      KSO3 (part 3) 
AB2C2DE3_hR9_155_b_c_c_a_e-001.{ABCDE}                               hR9       R32                 -             BBe2F2KO3               #155      KBe2BO3F2 (part 3) 
ABC3_hR5_160_a_a_b-001.{ABC}                                         hR5       R3m                 G0_{7}        BrKO3                   #160      KBrO3 (part 3) 
ABC3_hR5_160_a_a_b-002.{ABC}                                         hR5       R3m                 -             KNO3                    #160      gamma-Potassium Nitrate (KNO3) (part 3) 
AB3C2D_hR7_160_a_b_2a_a-001.{ABCD}                                   hR7       R3m                 S5_{7}        FeO3(OH)2Si             #160      Cronstedtite Fe(Fe,Si)[(OH)2,O]O3 (part 3) 
A3B7C_hR11_160_b_a2b_a-001.{ABC}                                     hR11      R3m                 -             Fe3O7P                  #160      Fe3PO7 (part 3) 
AB4C8_hR13_160_a_ab_2a2b-001.{ABC}                                   hR13      R3m                 -             GaMo4S8                 #160      Low-Temperature GaMo4S8 (part 3) 
A3B24C_hR28_160_b_2b3c_a-001.{ABC}                                   hR28      R3m                 -             I3S24Sb                 #160      SbI3S24 (part 3) 
A2BC4_hR42_161_2b_b_4b-001.{ABC}                                     hR42      R3c                 -             B2BaO4                  #161      alpha-BaB2O4 (Low-Temperature) (part 3) 
A2B6C_hP9_162_d_k_a-001.{ABC}                                        hP9       P-31m               I1_{3}        Cl2(H2O)6Sr             #162      I13 (SrCl2.(H2O)6) (obsolete) (part 3) 
A6BC2_hP9_162_k_a_d-001.{ABC}                                        hP9       P-31m               -             O6PbSb2                 #162      Rosiaite (PbSb2O6) (part 3) 
A6BC_hP16_163_i_b_c-001.{ABC}                                        hP16      P-31c               J1_{12}       F4NaSb                  #163      NaSbF4(OH)2 (part 3) 
ABC6D_hP18_163_d_b_i_c-001.ABCD/DACB/CABD/                           hP18      P-31c               -             AlCaF6Li                #163      Colquiriite (LiCaAlF6) (part 3) 
DCBA/BCAD/DBAC/CDAB/ADBC/BADC/CBDA/ACDB/BDCA
AB3_hP4_164_b_ad-001.{AB}                                            hP4       P-3m1               D0_{13}       AlCl3                   #164      D013 (AlCl3) (obsolete) (part 3) 
A2BC2_hP5_164_d_a_d-001.{ABC}                                        hP5       P-3m1               -             H2MgO2                  #164      Brucite [Mg(OH)2] (part 3) 
A2BC6_hP9_164_d_a_i-001.{ABC}                                        hP9       P-3m1               H6_{3}        K2Pt(SCN)6              #164      K2Pt(SCN)6 (part 3) 
A6B2C_hP9_164_i_d_a-001.{ABC}                                        hP9       P-3m1               J1_{6}        F6(NH4)2Si              #164      Bararite (Trigonal (NH4)2SiF6) (part 3) 
A6BC2_hP9_164_i_a_d-001.{ABC}                                        hP9       P-3m1               J1_{13}       F6GeK2                  #164      K2GeF6 (part 3) 
AB2C3_hP12_164_d_ae_i-001.{ABC}                                      hP12      P-3m1               -             HgPt2Se3                #164      Jacutingaite (Pt2HgSe3) (part 3) 
AB_hP12_164_c2d_c2d-001.{AB}                                         hP12      P-3m1               -             BiSe                    #164      Nevskite (BiSe) (part 3) 
A16B2C_hP19_164_2d2i_d_b-001.{ABC}                                   hP19      P-3m1               -             H16Li2Mg                #164      Predicted Li2MgH16 300GPa (part 3) 
AB_hR2_166_a_b-003.AB                                                hR2       R-3m                B22           KS                      #166      K(SH) (part 3) 
AB2_hR6_166_c_2c-002.{AB}                                            hR6       R-3m                C12           CaSi2                   #166      CaSi2 (part 3) 
AB4C_hR6_166_b_2c_a-001.{ABC}                                        hR6       R-3m                -             CaO4U                   #166      CaUO4 (part 3) 
A2BC4_hR7_166_c_a_2c-002.{ABC}                                       hR7       R-3m                -             Bi2MnTe4                #166      MnBi2Te4 (part 3) 
AB4C2_hR7_166_a_2c_c-001.{ABC}                                       hR7       R-3m                -             CaCu4P2                 #166      CaCu4P2 (part 3) 
ABC6_hR8_166_a_b_h-001.{ABC}                                         hR8       R-3m                -             KNO3                    #166      beta-Potassium Nitrate (KNO3) (part 3) 
A13B2_hR15_166_b2h_c-001.{AB}                                        hR15      R-3m                D1_{g}        B13C2                   #166      B13C2 ``B4C'' (part 3) 
AB4C2_hR28_166_2c_2c2h_abh-001.{ABC}                                 hR28      R-3m                -             CuS4Ti2                 #166      Rhombohedral CuTi2S4 (part 3) 
A5B21C24D12_hR62_166_a2c_ehi_fg2h_i-001.{ABCD}                       hR62      R-3m                S3_{4}(I)     Ca5(H2O)21O24Si12       #166      Chabazite (Ca1.4Sr0.3Al3.8Si8.3O24.13H2O) (part 3) 
A3B_hR8_167_e_b-001.{AB}                                             hR8       R-3c                D0_{12}       F3Fe                    #167      FeF3 (part 3) equivalent to 443,444 [HTQC]
A3BC6_hR20_167_e_b_f-001.{ABC}                                       hR20      R-3c                J2_{2}        Cl3Cr(H2O)6             #167      CrCl3(H2O)6 (part 3) 
A6BC3D_hR22_167_f_b_e_a-001.{ABCD}                                   hR22      R-3c                -             Cl6FeK3Na               #167      Rinneite (K3NaFeCl6) (part 3) 
A9B3C2_hR28_167_ef_e_c-001.{ABC}                                     hR28      R-3c                K7_{2}        Cl9Cs3Tl2               #167      Cs3Tl2Cl9 (part 3) 
A2BC4_hR42_167_f_ac_2f-001.{ABC}                                     hR42      R-3c                -             B2BaO4                  #167      beta-BaB2O4 (High-Temperature) (part 3) 
A25B21_hR92_167_b2e3f_e3f-001.{AB}                                   hR92      R-3c                -             Re25Zr21                #167      Zr21Re25 (part 3) equivalent to 481,482 [HTQC]
ABC3_hP10_173_b_a_c-001.{ABC}                                        hP10      P6_{3}              -             ILiO3                   #173      alpha-LiIO3 (part 3) 
ABC4D_hP14_173_a_b_bc_b-001.{ABCD}                                   hP14      P6_{3}              H1_{4}        KLiO4S                  #173      LiKSO4 (part 3) 
AB3C7D_hP24_173_a_c_b2c_b-001.{ABCD}                                 hP24      P6_{3}              -             CuLa3S7Si               #173      La3CuSiS7 (part 3) 
AB3C9D_hP28_173_a_c_3c_b-001.{ABCD}                                  hP28      P6_{3}              -             BLa3O9W                 #173      La3BWO9 (P6_{3}) (part 3) 
A3BCD3E15F3_hP52_173_c_b_b_c_5c_c-001.{ABCDEF}                       hP52      P6_{3}              S3_{3}(I)     Al3CCaNa3O15Si3         #173      Crancrinite (Na6Ca2Al6Si6O24(CO3)2) (part 3) 
A3B2_hP20_176_2h_ah-001.{AB}                                         hP20      P6_{3}/m            D8_{k}        S12Th7                  #176      Th7S12 (part 3) 
A9B3C2_hP28_176_hi_af_f-001.{ABC}                                    hP28      P6_{3}/m            K7_{1}        Cl9K3W2                 #176      K3W2Cl9 (part 3) 
A10B7_hP34_176_c3h_b2h-001.{AB}                                      hP34      P6_{3}/m            -             Rh20Si13                #176      Rh20Si13 (part 3) 
A9B2C9_hP40_176_hi_f_hi-001.{ABC}                                    hP40      P6_{3}/m            F4_{1}        C9Fe2O9                 #176      Fe2(CO)9 (part 3) 
A5BC12D3_hP42_176_fh_a_2hi_h-001.{ABCD}                              hP42      P6_{3}/m            H5_{7}        Ca5FO12P3               #176      Fluorapatite [Ca5F(PO4)3] (part 3) 
A2BCD2_hP18_180_f_c_b_i-001.{ABCD}                                   hP18      P6_{2}22            -             Hg2INaO2                #180      Hg2O2NaI (part 3) 
ABC3_hP10_182_c_b_g-001.{ABC}                                        hP10      P6_{3}22            E2_{3}        ILiO3                   #182      E23 (LiIO3) (obsolete) (part 3) 
A2BC6_hP18_182_f_b_gh-001.{ABC}                                      hP18      P6_{3}22            H2_{8}        Al2BaO4                 #182      BaAl2O4 (part 3) 
AB2_hP6_186_b_ab-001.{AB}                                            hP6       P6_{3}mc            C27           CdI2                    #186      C27 (CdI2) (emquestionable) (part 3) 
ABCD_hP8_186_b_b_a_a-001.{ABCD}                                      hP8       P6_{3}mc            E0_{3}        CdClHO                  #186      Cd(OH)Cl (part 3) 
A3B8C2_hP26_186_c_ab2c_2b-001.{ABC}                                  hP26      P6_{3}mc            -             Mo3O8Zn2                #186      Zn2Mo3O8 (part 3) 
A4BC7D_hP26_186_ac_b_a2c_b-001.{ABCD}                                hP26      P6_{3}mc            E9_{2}        Be4NaO7Sb               #186      Swedenborgite (NaBe4SbO7) (part 3) 
AB6CD7_hP30_186_b_d_a_b2c-001.{ABCD}                                 hP30      P6_{3}mc            H4_{18}       ClH6LiO7                #186      LiClO4.3H2O (part 3) 
A3B9CD9_hP44_186_c_3c_b_cd-001.{ABCD}                                hP44      P6_{3}mc            G2_{2}        Br3(H2O)9NdO9           #186      Nd(BrO3)3.9H2O (part 3) 
A3B3C2_hP16_187_jk_jk_ck-001.{ABC}                                   hP16      P-6m2               -             As3Cr3K2                #187      Cr-233 Quasi-One-Dimensional Superconductor (K2Cr3As3) (part 3) 
A7B_hP24_187_ai2j2kn_j-001.{AB}                                      hP24      P-6m2               -             Cs7O                    #187      Cs7O (part 3) 
AB3C_hP20_190_ac_i_f-001.{ABC}                                       hP20      P-62c               K1_{2}        CsO3S                   #190      CsSO3 (part 3) 
ABCD3_hP36_190_h_g_af_hi-001.{ABCD}                                  hP36      P-62c               G7_{1}        CCeFO3                  #190      Bastn"asite [CeF(CO3)] (part 3) 
A12B_hP13_191_cdei_a-001.{AB}                                        hP13      P6/mmm              D2_{a}        Be12Ti                  #191      TiBe12 (approximate) (part 3) 
A4B5_hP18_193_bg_dg-001.{AB}                                         hP18      P6_{3}/mcm          -             Ga4Ti5                  #193      Ti5Ga4 (part 3) equivalent to 359,360 [HTQC]
A3B_hP24_193_ack_g-001.{AB}                                          hP24      P6_{3}/mcm          D0_{6}        F3La                    #193      D06 (Tysonite, LaF3) (obsolete) (part 3) 
AB_hP4_194_a_c-003.{AB}                                              hP4       P6_{3}/mmc          C_{k}         LiZn2                   #194      LiZn2 (part 3) 
A2B3_hP10_194_f_bf-001.{AB}                                          hP10      P6_{3}/mmc          D5_{b}        Pt2Sn3                  #194      Pt2Sn3 (part 3) 
AB2C2_hP10_194_a_bc_f-001.{ABC}                                      hP10      P6_{3}/mmc          -             CoNa0.74O2              #194      Na0.74CoO2 (part 3) 
AB2C2_hP10_194_a_f_f-001.{ABC}                                       hP10      P6_{3}/mmc          -             EuIn2P2                 #194      EuIn2P2 (part 3) 
A10B_hP22_194_bhj_c-001.{AB}                                         hP22      P6_{3}/mmc          -             H10Hf                   #194      Proposed 300GPa HfH10 (part 3) 
AB3C2_hP24_194_f_k_bh-001.{ABC}                                      hP24      P6_{3}/mmc          -             CoGa3Lu2                #194      Lu2CoGa3 (part 3) 
A9B2C3_hP28_194_hk_f_bf-001.{ABC}                                    hP28      P6_{3}/mmc          -             Cl9Cr2Cs3               #194      Cs3Cr2Cl9 (part 3) 
A3B2C9D3E_hP36_194_g_f_hk_h_a-001.{ABCDE}                            hP36      P6_{3}/mmc          S3_{4}(II)    (H2O)Na2O9Si3Zr         #194      S34 (II) (Catapleiite, Na2Zr(SiO3)3.H2O) (obsolete) (part 3) 
A2B3_hP60_194_3fk_cdef2k-001.{AB}                                    hP60      P6_{3}/mmc          D5_{6}        Al2O3                   #194      beta-Alumina (Al2O3) (part 3) 
A12B19C_hP64_194_ab2fk_efh2k_d-001.{ABC}                             hP64      P6_{3}/mmc          -             Fe12O19Pb               #194      Magnetoplumbite (PbFe12O19) (part 3) 
ABC4D_cP28_198_a_a_ab_a-001.{ABCD}                                   cP28      P2_{1}3             S6_{5}        AlNaO4Si                #198      alpha-Carnegieite (NaAlSiO4) (part 3) 
AB2C4D_cP32_198_a_2a_ab_a-001.{ABCD}                                 cP32      P2_{1}3             S6_{6}        CaNa2O4Si               #198      Na2CaSiO4 (part 3) 
A2B4C_cP56_198_ab_2a2b_2a-001.{ABC}                                  cP56      P2_{1}3             -             Cu2O4Se                 #198      Cubic Cu2OSeO3 (part 3) 
AB2_cI36_199_b_c-001.{AB}                                            cI36      I2_{1}3             C26_{a}       NO2                     #199      C26a (NO2) (obsolete) (part 3) 
A3B11C3_cP68_201_be_efh_g-001.{ABC}                                  cP68      Pn-3                -             Bi3O11Ru3               #201      Bi3Ru3O11 (part 3) 
AB3C6D12_cF88_202_a_bc_e_h-001.{ABCD}                                cF88      Fm-3                J2_{4}        CoK3N6O12               #202      K3Co(NO2)6 (part 3) 
A4BC12_cI34_204_c_a_g-001.{ABC}                                      cI34      Im-3                -             Fe4LaP12                #204      LaFe4P12 (part 3) 
AB2_cI36_204_d_g-001.{AB}                                            cI36      Im-3                C26           NO2                     #204      NO2 (Modern) (part 3) 
A7BC12_cI40_204_bc_a_g-001.{ABC}                                     cI40      Im-3                -             Mn7NaO12                #204      NaMn7O12 (part 3) 
A6BC_cP32_205_d_b_a-001.{ABC}                                        cP32      Pa-3                -             F6NaSb                  #205      NaSbF6 (part 3) 
A2B6C_cP36_205_c_d_a-001.{ABC}                                       cP36      Pa-3                G2_{1}        N2O6Pb                  #205      Pb(NO3)2 (part 3) 
A4B_cP40_205_cd_c-001.{AB}                                           cP40      Pa-3                D1_{1}        I4Sn                    #205      SnI4 (part 3) 
A7B2C_cP40_205_bd_c_a-001.{ABC}                                      cP40      Pa-3                K6_{1}        P2O7Zr                  #205      ZrP2O7 High-Temperature (part 3) 
A2B6C6D_cP60_205_c_d_d_a-001.{ABCD}                                  cP60      Pa-3                J1_{10}       Br2(H2O)6O6Zn           #205      Zn(BrO3)2.6H2O (J110) (part 3) 
A2B6CD6_cP60_205_c_d_a_d-001.{ABCD}                                  cP60      Pa-3                H6_{4}        N2(NH3)6NiO6            #205      H64 [Ni(NO3)2(NH3)6] (obsolete) (part 3) 
A2BC4_cP84_205_d_ac_2d-001.{ABC}                                     cP84      Pa-3                -             B2CaO4                  #205      CaB2O4 (IV) (part 3) 
AB12CD8E2_cP96_205_a_2d_b_cd_c-001.{ABCDE}                           cP96      Pa-3                -             Cr(H2O)12NaO8S2         #205      NaCr(SO4)2.12H2O Alum (part 3) 
AB24CD20E2_cP192_205_a_4d_b_c3d_c-001.{ABCDE}                        cP192     Pa-3                H4_{15}       AlH24NaO20S2            #205      gamma-Alum [AlNa(SO4)2.12H2O] (part 3) 
AB24CD28E2_cP224_205_a_4d_b_2c4d_c-001.{ABCDE}                       cP224     Pa-3                H4_{13}       AlH24KO20S2             #205      alpha-Alum [KAl(SO4)2.12H2O] (part 3) 
AB2C36D2E20F2_cP252_205_a_c_6d_c_c3d_c-001.{ABCDEF}                  cP252     Pa-3                H4_{14}       AlC2H36N2O20S2          #205      beta-Alum [Al(NH3CH3)2(SO4)2.12H2O] (part 3) 
A2B3_cP60_212_bcd_ace-001.{AB}                                       cP60      P4_{3}32            D5_{7}        Fe2O3                   #212      Maghemite (gamma-Fe2O3) (part 3) 
A3B2_cP20_213_d_c-001.{AB}                                           cP20      P4_{1}32            -             Mg3Ru2                  #213      Mg3Ru2 (part 3) equivalent to 593,594 [HTQC]
A2BC3_cP24_213_c_a_d-001.{ABC}                                       cP24      P4_{1}32            -             Al2CMo3                 #213      Al2Mo3C (part 3) 
ABC4_cF24_216_b_a_e-001.{ABC}                                        cF24      F-43m               H0_{5}        ClKO4                   #216      High-Temperature Cubic KClO4 (part 3) 
AB_cF40_216_de_ce-001.AB                                             cF40      F-43m               -             AlN                     #216      AlN (cF40) (part 3) 
AB4C8_cF52_216_a_e_2e-001.{ABC}                                      cF52      F-43m               -             GaMo4S8                 #216      GaMo4S8 (part 3) 
A13BC18D20E5_cF228_216_dh_b_fh_2eh_ce-001.{ABCDE}                    cF228     F-43m               S0_{8}        Al13ClF18O20Si5         #216      Zunyite [Al13(OH,F)18Si5O20Cl] (part 3) 
A16B40C12D6E5_cF316_216_eh_e2g2h_h_f_be-001.{ABCDE}                  cF316     F-43m               -             F16O40Ti12Y6Zn5         #216      Murataite [(Y,Na)6(Zn,Fe)5Ti12O29(O,F)10F4] (part 3) 
A45B11_cF448_216_bd4efg5h_ac2eh-001.{AB}                             cF448     F-43m               -             Cd45Sm11                #216      Sm11Cd45 (part 3) 
AB_cI16_217_c_c-001.AB                                               cI16      I-43m               -             AlN                     #217      AlN (cI16) (part 3) 
A4B24C13_cI82_217_c_deg_ag-001.{ABC}                                 cI82      I-43m               -             As4Cu12S13              #217      Tennantite (Cu12As4S13) (part 3) 
A3BC4D12E3_cP46_218_d_a_e_i_c-001.{ABCDE}/2 (60)                     cP46      P-43n               S6_{2}        Al3ClNa4O12Si3          #218      Sodalite [Na4(AlSiO4)3Cl] (part 3) 
A3B4C4D4E16F4G3_cP76_218_c_e_e_e_ei_e_d-001.{ABCDEFG}/2 (2520)       cP76      P-43n               S6_{9}        Al3Ca4K4Na4O16S4Si3     #218      Hauyne [(Na0.5Ca0.3K0.2)8(Al6Si6O24)(SO4)1.5] (part 3) 
AB_cI24_220_a_b.AB                                                   cI24      I-43d               -             AlN                     #220      AlN (cI24) (part 3) 
A4B12C3_cI76_220_c_e_a-001.{ABC}                                     cI76      I-43d               S1_{5}        Bi4O12Si3               #220      Eulytine (Bi4(SiO4)3) (part 3) 
A7B12C19_cI152_220_bc_2d_ace-001.{ABC}                               cI152     I-43d               K7_{4}        Al14Ca12O33             #220      Mayenite (12CaO.7Al2O3, C12A7) (part 3) 
AB9C3_cI208_220_c_3e_e-001.{ABC}                                     cI208     I-43d               G5_{2}        AlO9P3                  #220      Al(PO3)3 (part 3) 
AB12C_cP14_221_a_h_b-001.ABC/CAB/BCA                                 cP14      Pm-3m               -             CaH12Y                  #221      Predicted High-Pressure YCaH12 (part 3) 
A3B40CD12_cP112_224_d_e3k_a_k-001.{ABCD}                             cP112     Pn-3m               -             (H3O)3O40PW12           #224      H3PW12O40.3H2O (part 3) 
A5B40CD12_cP116_224_cd_e3k_a_k-001.{ABCD}                            cP116     Pn-3m               H4_{16}       (H2.6O)5O40PW12         #224      12-phosphotungstic acid [H3PW12O40.5H2O] (part 3) 
A27B52CD12_cP184_224_dl_eh3k_a_k-001.{ABCD}                          cP184     Pn-3m               -             H15O46PW12              #224      Dodecatungstophosphoric Acid Hexahydrate [H3PW12O40.6H2O] (part 3) 
A2BC6D_cF40_225_c_a_e_b-001.ABCD/BDCA/CABD/                          cF40      Fm-3m               -             Ba2MnO6W                #225      Double Perovskite (Ba2MnWO6) (part 3) equivalent to Q0001 [HTQC]
ADBC/BCAD/CDAB/DBAC/DCBA/BADC/ACDB/CBDA/DACB
A10B_cF44_225_cf_b-001.{AB}                                          cF44      Fm-3m               -             H10La                   #225      LaH10 High-Tc Superconductor (part 3) 
A9B8_cF68_225_af_ce-001.{AB}                                         cF68      Fm-3m               D8_{9}        Co9S8                   #225      Co9S8 (part 3) 
ABC6D8E2_cF72_225_b_a_e_f_c-001.{ABCDE}                              cF72      Fm-3m               H5_{8}        ClFNa6O8S2              #225      Sulphohalite [Na6ClF(SO4)2] (part 3) 
A6B9CD2E6_cF96_225_e_bf_a_c_e-001.{ABCDE}                            cF96      Fm-3m               J2_{5}        C12Cu3Fe2(H2O)xN12      #225      Cu3[Fe(CN)6]2. xH2O (x approx 3) (part 3) 
AB30C16D3_cF200_225_a_ej_2f_bc-001.{ABCD}                            cF200     Fm-3m               J2_{1}        AlF6H12N3               #225      (NH4)3AlF6 (part 3) 
AB_cF32_227_c_d.AB                                                   cF32      Fd-3m               L1_{3}(I)     CuPt                    #227      Cubic CuPt (part 3) 
A3B2_cF80_227_f_e-002.{AB}                                           cF80      Fd-3m               D6_{1}        O3Sb2                   #227      Senarmontite (Sb2O3) (part 3) 
A2B_cF96_227_abf_cd-001.{AB}                                         cF96      Fd-3m               D6_{2}        O2Sb                    #227      D62 (Sb2O4) (obsolete) (part 3) 
A11B4_cF120_227_acdf_e-001.{AB}                                      cF120     Fd-3m               -             Ga2O3                   #227      gamma-Ga2O3 (part 3) 
A16B2C_cF152_227_eg_d_a-001.{ABC}                                    cF152     Fd-3m               -             H16Li2Mg                #227      Predicted Li2MgH16 High-Temperature Superconductor (250GPa) (part 3) 
A18B2C3_cF184_227_fg_d_ac-001.{ABC}                                  cF184     Fd-3m               -             Al18Cr2Mg3              #227      Mg3Cr2Al18 (part 3) 
A22B_cF184_227_cdfg_a-001.{AB}                                       cF184     Fd-3m               -             Zn22Zr                  #227      Zn22Zr (part 3) equivalent to 425,426 [HTQC]
A2BCD3E6_cF208_227_e_c_d_f_g-001.{ABCDE}                             cF208     Fd-3m               G7_{3}        C2ClMgNa3O6             #227      G73 [Northupite, Na3MgCl(CO3)2] (obsolete) (part 3) 
A4B2C6D16E_cF232_227_e_d_f_eg_a-001.{ABCDE}                          cF232     Fd-3m               H5_{6}        C4Mg2Na6O16S            #227      H56 [Tychite, Na6Mg2SO4(CO3)4)] (obsolete) (part 3) 
A29B40CD12_cF656_227_ae2fg_e3g_b_g-001.{ABCD}                        cF656     Fd-3m               H4_{21}       (H2O)29O40PW12          #227      H3PW12O40.29H2O (part 3) 
A21B_cI44_229_bdh_a-001.{AB}                                         cI44      Im-3m               B23           AgI                     #229      alpha-AgI (part 3) 
A2B3C12D12_cI232_230_a_c_h_h-001.{ABCD}                              cI232     Ia-3d               J2_{3}        Al2Ca3H12O12            #230      Ca3Al2(OH)12 (part 3) 

*****************************************************************************************************************************************************
ANRL - metal carbo-nitrides - QUATERNARIES

ANRL Label                                                           Pearson   SG symbol           Struk         Prototype               SG #      notes
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A2BCD4_mP16_11_2e_e_e_4e-001.{ABCD}                                  mP16      P2_{1}/m            -             C2La1Li1N4              #11       quaternary metal-carbo-nitride prototype; C2La1Li1N4, ICSD #420440
AB2C2D_mC24_15_a_f_f_e-001.{ABCD}                                    mC24      C2/c                -             Ag1C2N2Na1              #15       quaternary metal-carbo-nitride prototype; Ag1C2N2Na1, ICSD #65697
A6BC6D2_hP15_147_g_a_g_d-001.{ABCD}                                  hP15      P-3                 -             C6Fe1N6Pb2              #147      quaternary metal-carbo-nitride prototype; C6Fe1N6Pb2, ICSD #51493
A6BC6D2_hP15_147_g_a_g_d-002.{ABCD}                                  hP15      P-3                 -             C6Fe1N6Sn2              #147      quaternary metal-carbo-nitride prototype; C6Fe1N6Sn2, ICSD #51494
A6B2C6D_hP15_147_g_d_g_a-001.{ABCD}                                  hP15      P-3                 -             C6Cd2N6Os1              #147      quaternary metal-carbo-nitride prototype; C6Cd2N6Os1, ICSD #417821
A6B2C6D_hP15_147_g_d_g_a-002.{ABCD}                                  hP15      P-3                 -             C6Cd2N6Ru1              #147      quaternary metal-carbo-nitride prototype; C6Cd2N6Ru1, ICSD #417822
A6B2C6D_hP15_147_g_d_g_a-003.{ABCD}                                  hP15      P-3                 -             C6Mn2N6Ru1              #147      quaternary metal-carbo-nitride prototype; C6Mn2N6Ru1, ICSD #417823
A6B2CD6_hP15_147_g_d_a_g-001.{ABCD}                                  hP15      P-3                 -             C6Cd2Fe1N6              #147      quaternary metal-carbo-nitride prototype; C6Cd2Fe1N6, ICSD #417826
A3B6CD6_hP16_162_g_k_a_k-001.{ABCD}                                  hP16      P-31m               -             Ag3C6Co1N6              #162      quaternary metal-carbo-nitride prototype; Ag3C6Co1N6, ICSD #16959
A3B6CD6_hP16_162_g_k_a_k-002.{ABCD}                                  hP16      P-31m               -             Ag3C6Fe1N6              #162      quaternary metal-carbo-nitride prototype; Ag3C6Fe1N6, ICSD #173552
A6B2CD6_cF60_216_f_bc_a_f-001.{ABCD}                                 cF60      F-43m               -             C6Cu2Fe1N6              #216      quaternary metal-carbo-nitride prototype; C6Cu2Fe1N6, ICSD #28653
A6BC6D_cF56_225_e_a_e_b-001.{ABCD}                                   cF56      Fm-3m               -             C6Cd1N6Pd1              #225      quaternary metal-carbo-nitride prototype; C6Cd1N6Pd1, ICSD #6093
A3BCD3_cF64_225_e_c_ab_e-001.{ABCD}                                  cF64      Fm-3m               -             C6Cs2Mn2N6              #225      quaternary metal-carbo-nitride prototype; C6Cs2Mn2N6, ICSD #248046

****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
