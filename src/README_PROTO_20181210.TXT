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
1001.AB/BA     cF52     Fm-3m                    B12Y1           #225    ICSD_23860    binary_boride_ICSD_frequency=16
1002.AB/BA     tP20     P4/mbm                   B4La1           #127    ICSD_2360     binary_boride_ICSD_frequency=14
1003.AB/BA     hP3      P6/mmm                   B2Nb1           #191    ICSD_30328    binary_boride_ICSD_frequency=10
1004.AB/BA     hP3      P6/mmm                   B2Pu1           #191    ICSD_43660    binary_boride_ICSD_frequency=10
1005.AB/BA     tI12     I4/mcm                   B1Mo2           #140    ICSD_24278    binary_boride_ICSD_frequency=8
1006.AB/BA     oS8      Cmcm                     B1Ni1           #63     ICSD_26937    binary_boride_ICSD_frequency=8
1007.AB/BA     hP3      P6/mmm                   B2Mg1           #191    ICSD_26675    binary_boride_ICSD_frequency=6
1008.AB        cF8      F-43m                    As1B1           #216    ICSD_43871    binary_boride_ICSD_frequency=5
1009.AB/BA     oI14     Immm                     B4Cr3           #71     ICSD_33704    binary_boride_ICSD_frequency=5
1010.AB        cF8      Fm-3m                    B1Pu1           #225    ICSD_43661    binary_boride_ICSD_frequency=4
1011.AB/BA     tP10     P4/mbm                   B2V3            #127    ICSD_88317    binary_boride_ICSD_frequency=4
1012.AB/BA     oP8      Pnma                     B1Ti1           #62     ICSD_24701    binary_boride_ICSD_frequency=3
1013.AB/BA     hP20     P6_{3}mc                 B3Ru7           #186    ICSD_44343    binary_boride_ICSD_frequency=3
1014.AB/BA     oS20     Cmcm                     B3V2            #63     ICSD_79258    binary_boride_ICSD_frequency=3
1015.AB/BA     cF116    Fm-3m                    B6Fe23          #225    ICSD_54786    binary_boride_ICSD_frequency=3
1016.AB/BA     hP3      P6/mmm                   Ag1B2           #191    ICSD_43821    binary_boride_ICSD_frequency=2
1017.AB/BA     oP8      Pnma                     B1Co1           #62     ICSD_612863   binary_boride_ICSD_frequency=2
1018.AB/BA     tI12     I4/mcm                   B1Co2           #140    ICSD_42531    binary_boride_ICSD_frequency=2
1019.AB/BA     oP16     Pnma                     B1Co3           #62     ICSD_44339    binary_boride_ICSD_frequency=2
1020.AB/BA     tI16     I4_{1}/amd               B1Mo1           #141    ICSD_24280    binary_boride_ICSD_frequency=2
1021.AB/BA     oS16     Cmcm                     B1Re3           #63     ICSD_43662    binary_boride_ICSD_frequency=2
1022.AB        hP2      P-6m2                    B1Tc1           #187    ICSD_168895   binary_boride_ICSD_frequency=2
1023.AB/BA     hP6      P6_{3}/mmc               B2Re1           #194    ICSD_23871    binary_boride_ICSD_frequency=2
1024.AB/BA     oP6      Pmmn                     B2Ru1           #59     ICSD_31871    binary_boride_ICSD_frequency=2
1025.AB/BA     hP10     P6_{3}/mmc               B3Ru2           #194    ICSD_23715    binary_boride_ICSD_frequency=2
1026.AB/BA     oP10     Pnnm                     B4Cr1           #58     ICSD_186851   binary_boride_ICSD_frequency=2
1027.AB/BA     oI10     C2/m                     B4Mn1           #12     ICSD_15079    binary_boride_ICSD_frequency=2
1028.AB/BA     hR7      R-3m                     B5W2            #166    ICSD_20326    binary_boride_ICSD_frequency=2
1029.AB/BA     cF48     F-43m                    B11Li1          #216    ICSD_164841   binary_boride_ICSD_frequency=1
1030.AB/BA     hR12     R3m                      B11Li1          #160    ICSD_164842   binary_boride_ICSD_frequency=1
1031.AB/BA     cF48     F-43m                    B11Li1          #216    ICSD_164843   binary_boride_ICSD_frequency=1
1032.AB/BA     tI24     I-4m2                    B11Li1          #119    ICSD_164844   binary_boride_ICSD_frequency=1
1033.AB/BA     tI24     I-4m2                    B11Li1          #119    ICSD_164845   binary_boride_ICSD_frequency=1
1034.AB/BA     tI26     I4/mmm                   B12Sc1          #139    ICSD_615424   binary_boride_ICSD_frequency=1
1035.AB/BA     hR15     R-3m                     B12Si3          #166    ICSD_615435   binary_boride_ICSD_frequency=1
1036.AB/BA     cF12     Fm-3m                    B1Be2           #225    ICSD_20384    binary_boride_ICSD_frequency=1
1037.AB/BA     oP16     Pnma                     B1Co3           #62     ICSD_603543   binary_boride_ICSD_frequency=1
1038.AB/BA     oP16     Pnma                     B1Co3           #62     ICSD_612862   binary_boride_ICSD_frequency=1
1039.AB/BA     oF48     Fddd                     B1Cr2           #70     ICSD_188474   binary_boride_ICSD_frequency=1
1040.AB/BA     tP3      I4/mcm                   B1Fe2           #140    ICSD_160789   binary_boride_ICSD_frequency=1
1041.AB/BA     tI12     I-42m                    B1Fe2           #121    ICSD_16809    binary_boride_ICSD_frequency=1
1042.AB/BA     hP8      P6_{3}22                 B1Fe3           #182    ICSD_184958   binary_boride_ICSD_frequency=1
1043.AB/BA     oP8      Pnma                     B1Li1           #62     ICSD_153291   binary_boride_ICSD_frequency=1
1044.AB/BA     oP16     Pnma                     B1Ni3           #62     ICSD_24306    binary_boride_ICSD_frequency=1
1045.AB/BA     oP6      Pnnm                     B1Pd2           #58     ICSD_10487    binary_boride_ICSD_frequency=1
1046.AB/BA     oP16     Pnma                     B1Pd3           #62     ICSD_43514    binary_boride_ICSD_frequency=1
1047.AB/BA     hP4      P6_{3}/mmc               B1Pt1           #194    ICSD_24363    binary_boride_ICSD_frequency=1
1048.AB/BA     hP4      P6_{3}/mmc               B1Pt1           #194    ICSD_615210   binary_boride_ICSD_frequency=1
1049.AB/BA     hP6      P6_{3}/mmc               B1Pt2           #194    ICSD_615207   binary_boride_ICSD_frequency=1
1050.AB/BA     oS8      Cmcm                     B1Rh1           #63     ICSD_150732   binary_boride_ICSD_frequency=1
1051.AB/BA     oP12     Pnma                     B1Rh2           #62     ICSD_24699    binary_boride_ICSD_frequency=1
1052.AB/BA     hP4      P6_{3}/mmc               B1Tc1           #194    ICSD_168897   binary_boride_ICSD_frequency=1
1053.AB/BA     hP4      P6_{3}/mmc               B1Tc1           #194    ICSD_168898   binary_boride_ICSD_frequency=1
1054.AB        cP2      Pm-3m                    B1Tc1           #221    ICSD_168900   binary_boride_ICSD_frequency=1
1055.AB/BA     tI12     I4/mcm                   B1Ti2           #140    ICSD_189385   binary_boride_ICSD_frequency=1
1056.AB/BA     oS46     Cmmm                     B20Na3          #65     ICSD_407159   binary_boride_ICSD_frequency=1
1057.AB/BA     oP12     Pnma                     B2Fe1           #62     ICSD_425310   binary_boride_ICSD_frequency=1
1058.AB/BA     hR4      R-3m                     B2Li2           #166    Kolmogorov-Curtarolo LiB    binary_boride_ICSD_frequency=1
1059.AB/BA     hP8      P6_{3}/mmc               B2Li2           #194    Kolmogorov-Curtarolo LiB    binary_boride_ICSD_frequency=1
1060.AB/BA     mS28     C2/c                     B2Pd5           #15     ICSD_43513    binary_boride_ICSD_frequency=1
1061.AB/BA     cF24     Fd-3m                    B2Sc1           #227    ICSD_106903   binary_boride_ICSD_frequency=1
1062.AB/BA     oS6      P6/mmm                   B2Tm1           #191    ICSD_44598    binary_boride_ICSD_frequency=1
1063.AB/BA     hP12     P6_{3}/mmc               B2W1            #194    ICSD_23716    binary_boride_ICSD_frequency=1
1064.AB/BA     tI32     I4/mcm                   B3Cr5           #140    ICSD_41932    binary_boride_ICSD_frequency=1
1065.AB/BA     hR8      R-3m                     B3Mo1           #166    ICSD_167734   binary_boride_ICSD_frequency=1
1066.AB/BA     oP28     Pnma                     B3Ni4           #62     ICSD_24307    binary_boride_ICSD_frequency=1
1067.AB/BA     mS28     C2/c                     B3Ni4           #15     ICSD_24308    binary_boride_ICSD_frequency=1
1068.AB/BA     hP8      P6_{3}/mmc               B3Re1           #194    ICSD_24361    binary_boride_ICSD_frequency=1
1069.AB/BA     oS20     P6_{3}/mmc               B3Ru2           #194    ICSD_108082   binary_boride_ICSD_frequency=1
1070.AB/BA     oI64     Imma                     B3Si1           #74     ICSD_412621   binary_boride_ICSD_frequency=1
1071.AB/BA     oI14     Immm                     B4Mn3           #71     ICSD_44446    binary_boride_ICSD_frequency=1
1072.AB/BA     tP16     P4_{2}/mmc               B4Mn4           #131    ICSD_74985    binary_boride_ICSD_frequency=1
1073.AB/BA     hP10     P6_{3}/mmc               B4Mo1           #194    ICSD_167735   binary_boride_ICSD_frequency=1
1074.AB/BA     hP20     P6_{3}/mmc               B4Mo1           #194    ICSD_182095   binary_boride_ICSD_frequency=1
1075.AB/BA     hR6      R-3m                     B4Mo2           #166    ICSD_39554    binary_boride_ICSD_frequency=1
1076.AB/BA     hP18     P6_{3}/mmc               B4Rh5           #194    ICSD_86395    binary_boride_ICSD_frequency=1
1077.AB/BA     hP20     P6_{3}/mmc               B4W1            #194    ICSD_43193    binary_boride_ICSD_frequency=1
1078.AB/BA     tP20     P4/mbm                   B4W1            #127    ICSD_615683   binary_boride_ICSD_frequency=1
1079.AB/BA     mP28     P2_{1}/c                 B5Gd2           #14     ICSD_62067    binary_boride_ICSD_frequency=1
1080.AB/BA     hP14     P6_{3}/mmc               B5Mo2           #194    ICSD_167733   binary_boride_ICSD_frequency=1
1081.AB/BA     mS14     C2/m                     B5W2            #12     ICSD_20326    binary_boride_ICSD_frequency=1
1082.AB/BA     hP14     P6_{3}/mmc               B5W2            #194    ICSD_615696   binary_boride_ICSD_frequency=1
1083.AB/BA     oS22     Cmmm                     B6Ta5           #65     ICSD_68538    binary_boride_ICSD_frequency=1
1084.AB/BA     oP38     Pbam                     B8Ru11          #55     ICSD_43663    binary_boride_ICSD_frequency=1
1085.AB/BA     hP20     P6_{3}cm                 B9Li1           #185    ICSD_164625   binary_boride_ICSD_frequency=1


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
1200.AB        cF8      Fm-3m                    C1Ti1           #225    ICSD_1546     binary_carbide_ICSD_frequency=31
1201.AB        cF8      F-43m                    C1Si1           #216    ICSD_24217    binary_carbide_ICSD_frequency=20
1202.AB/BA     tI6      I4/mmm                   C2U1            #139    ICSD_14330    binary_carbide_ICSD_frequency=19
1203.AB/BA     cI40     I-43d                    C3Pr2           #220    ICSD_2444     binary_carbide_ICSD_frequency=17
1204.AB        hP2      P-6m2                    C1W1            #187    ICSD_15406    binary_carbide_ICSD_frequency=12
1205.AB        cP2      Pm-3m                    C1Sc1           #221    ICSD_181048   binary_carbide_ICSD_frequency=11
1206.AB/BA     hR3      R-3m                     C1Y2            #166    ICSD_22283    binary_carbide_ICSD_frequency=4
1207.AB/BA     tP70     P4/mnc                   C4Sc3           #128    ICSD_71145    binary_carbide_ICSD_frequency=4
1208.AB/BA     tP68     P-42_{1}c                C19Er15         #114    ICSD_76817    binary_carbide_ICSD_frequency=3
1209.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_16593    binary_carbide_ICSD_frequency=3
1210.AB/BA     hP3      P-3m1                    C1Nb2           #164    ICSD_43670    binary_carbide_ICSD_frequency=3
1211.AB/BA     mS12     C2/c                     C2Ca1           #15     ICSD_54184    binary_carbide_ICSD_frequency=3
1212.AB/BA     cP12     Pa-3                     C2Th1           #205    ICSD_77470    binary_carbide_ICSD_frequency=3
1213.AB/BA     oP18     Pbam                     C5Y4            #55     ICSD_84869    binary_carbide_ICSD_frequency=3
1214.AB/BA     cI132    Im-3                     C60K6           #204    ICSD_66879    binary_carbide_ICSD_frequency=3
1215.AB/BA     cF116    Fm-3m                    C6Cr23          #225    ICSD_62667    binary_carbide_ICSD_frequency=3
1216.AB/BA     mP22     P2_{1}/c                 C7Ho4           #14     ICSD_83381    binary_carbide_ICSD_frequency=3
1217.AB/BA     cF12     Fm-3m                    Be2C1           #225    ICSD_41567    binary_carbide_ICSD_frequency=2
1218.AB/BA     oP6      Pnnm                     C1Fe2           #58     ICSD_76826    binary_carbide_ICSD_frequency=2
1219.AB/BA     hP4      P6_{3}/mmc               C1Os1           #194    ICSD_168278   binary_carbide_ICSD_frequency=2
1220.AB/BA     oP12     Pbcn                     C1V2            #60     ICSD_9965     binary_carbide_ICSD_frequency=2
1221.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_43194    binary_carbide_ICSD_frequency=2
1222.AB/BA     hP3      P-6m2                    C2Os1           #187    ICSD_168279   binary_carbide_ICSD_frequency=2
1223.AB/BA     hP6      P6_{3}/mmc               C2Os1           #194    ICSD_168280   binary_carbide_ICSD_frequency=2
1224.AB/BA     oP16     Pnma                     C2Rb2           #62     ICSD_51529    binary_carbide_ICSD_frequency=2
1225.AB/BA     cI132    Im-3                     C60Sr6          #204    ICSD_75356    binary_carbide_ICSD_frequency=2
1226.AB/BA     hR7      R3m                      Al4C3           #160    ICSD_14397    binary_carbide_ICSD_frequency=1
1227.AB/BA     cP7      P-43m                    As4C3           #215    ICSD_163831   binary_carbide_ICSD_frequency=1
1228.AB/BA     hR15     R-3m                     B13C2           #166    ICSD_446      binary_carbide_ICSD_frequency=1
1229.AB/BA     tI24     I-4m2                    B1C5            #119    ICSD_166553   binary_carbide_ICSD_frequency=1
1230.AB/BA     oI12     Imm2                     B1C5            #44     ICSD_166554   binary_carbide_ICSD_frequency=1
1231.AB/BA     hP6      P3m1                     B1C5            #156    ICSD_166555   binary_carbide_ICSD_frequency=1
1232.AB/BA     aP12     P-1                      B1C5            #2      ICSD_166556   binary_carbide_ICSD_frequency=1
1233.AB/BA     tP8      P-4m2                    B1C7            #115    ICSD_181952   binary_carbide_ICSD_frequency=1
1234.AB/BA     hP8      P3m1                     B1C7            #156    ICSD_181953   binary_carbide_ICSD_frequency=1
1235.AB/BA     cP8      P-43m                    B1C7            #215    ICSD_181954   binary_carbide_ICSD_frequency=1
1236.AB/BA     oP8      Pmm2                     B1C7            #25     ICSD_181955   binary_carbide_ICSD_frequency=1
1237.AB/BA     hR8      R3m                      B1C7            #160    ICSD_181956   binary_carbide_ICSD_frequency=1
1238.AB/BA     hR15     R-3m                     B4C1            #166    ICSD_29093    binary_carbide_ICSD_frequency=1
1239.AB/BA     oP24     Pnma                     Ba1C2           #62     ICSD_168408   binary_carbide_ICSD_frequency=1
1240.AB/BA     hR3      R-3m                     Ba1C2           #166    ICSD_186576   binary_carbide_ICSD_frequency=1
1241.AB/BA     cP126    Pm-3n                    Ba3C60          #223    ICSD_70063    binary_carbide_ICSD_frequency=1
1242.AB/BA     cI132    Im-3                     C10Cs1          #204    ICSD_657337   binary_carbide_ICSD_frequency=1
1243.AB/BA     tP68     P-42_{1}c                C19Sc15         #114    ICSD_42631    binary_carbide_ICSD_frequency=1
1244.AB/BA     oP16     Pnma                     C1Cr3           #62     ICSD_603557   binary_carbide_ICSD_frequency=1
1245.AB/BA     hP3      P6/mmm                   C1Fe2           #191    ICSD_162103   binary_carbide_ICSD_frequency=1
1246.AB/BA     oP12     Pnma                     C1Fe2           #62     ICSD_187138   binary_carbide_ICSD_frequency=1
1247.AB/BA     tI32     I-4                      C1Fe3           #82     ICSD_187142   binary_carbide_ICSD_frequency=1
1248.AB/BA     hP8      P6_{3}22                 C1Fe3           #182    ICSD_42542    binary_carbide_ICSD_frequency=1
1249.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_43522    binary_carbide_ICSD_frequency=1
1250.AB/BA     mP10     P2_{1}/m                 C1Fe4           #11     ICSD_187143   binary_carbide_ICSD_frequency=1
1251.AB/BA     tI10     I4/m                     C1Fe4           #87     ICSD_187144   binary_carbide_ICSD_frequency=1
1252.AB/BA     cP5      P-43m                    C1Fe4           #215    ICSD_44729    binary_carbide_ICSD_frequency=1
1253.AB/BA     oP16     Pnma                     C1Mn3           #62     ICSD_603556   binary_carbide_ICSD_frequency=1
1254.AB/BA     cF112    Fd-3m                    C1Mn6           #227    ICSD_187037   binary_carbide_ICSD_frequency=1
1255.AB/BA     hP12     P6_{3}/mmc               C1Mo1           #194    ICSD_44987    binary_carbide_ICSD_frequency=1
1256.AB/BA     hP8      P6_{3}/mmc               C1Mo1           #194    ICSD_618298   binary_carbide_ICSD_frequency=1
1257.AB/BA     hP3      P-6m2                    C1Mo2           #187    ICSD_43669    binary_carbide_ICSD_frequency=1
1258.AB        hP4      P6_{3}mc                 C1Nb1           #186    ICSD_189092   binary_carbide_ICSD_frequency=1
1259.AB/BA     hP4      P6_{3}/mmc               C1Nb1           #194    ICSD_189093   binary_carbide_ICSD_frequency=1
1260.AB/BA     oP12     Pnma                     C1Nb2           #62     ICSD_31973    binary_carbide_ICSD_frequency=1
1261.AB/BA     oP12     Pmma                     C1Nb2           #51     ICSD_653912   binary_carbide_ICSD_frequency=1
1262.AB/BA     hR8      R-3c                     C1Ni3           #167    ICSD_17005    binary_carbide_ICSD_frequency=1
1263.AB/BA     cP8      P2_{1}3                  C1Os1           #198    ICSD_168277   binary_carbide_ICSD_frequency=1
1264.AB/BA     hP8      P6_{3}/mmc               C1Re1           #194    ICSD_618702   binary_carbide_ICSD_frequency=1
1265.AB/BA     hP6      P6_{3}/mmc               C1Re2           #194    ICSD_180376   binary_carbide_ICSD_frequency=1
1266.AB/BA     hR4      R-3m                     C1Ru1           #166    ICSD_188285   binary_carbide_ICSD_frequency=1
1267.AB/BA     hP3      P-3m1                    C1Sc2           #164    ICSD_280743   binary_carbide_ICSD_frequency=1
1268.AB        hP10	P3m1                     C1Si1           #156    ICSD_107204   binary_carbide_ICSD_frequency=1
1269.AB        hP12 	P6_{3}mc                 C1Si1           #186    ICSD_15325    binary_carbide_ICSD_frequency=1
1270.AB        hR16    	R3m                      C1Si1           #160    ICSD_18136    binary_carbide_ICSD_frequency=1
1271.AB        hR10    	R3m                      C1Si1           #160    ICSD_24168    binary_carbide_ICSD_frequency=1
1272.AB        hP8      P6_{3}mc                 C1Si1           #186    ICSD_24170    binary_carbide_ICSD_frequency=1
1273.AB        hP4      P6_{3}mc                 C1Si1           #186    ICSD_24261    binary_carbide_ICSD_frequency=1
1274.AB        hP16     P6_{3}mc                 C1Si1           #186    ICSD_24630    binary_carbide_ICSD_frequency=1
1275.AB        hR18     R3m                      C1Si1           #160    ICSD_24631    binary_carbide_ICSD_frequency=1
1276.AB        hP38     P3m1                     C1Si1           #156    ICSD_27635    binary_carbide_ICSD_frequency=1
1277.AB        hP28     P3m1                     C1Si1           #156    ICSD_28303    binary_carbide_ICSD_frequency=1
1278.AB        hP20     P3m1                     C1Si1           #156    ICSD_38150    binary_carbide_ICSD_frequency=1
1279.AB        hR14     R3m                      C1Si1           #160    ICSD_42513    binary_carbide_ICSD_frequency=1
1280.AB        hP32     P6_{3}mc                 C1Si1           #186    ICSD_42858    binary_carbide_ICSD_frequency=1
1281.AB        hP36     P3m1                     C1Si1           #156    ICSD_42859    binary_carbide_ICSD_frequency=1
1282.AB        hP42     P3m1                     C1Si1           #156    ICSD_43827    binary_carbide_ICSD_frequency=1
1283.AB/BA     cF48     Fd-3m                    C1Ti2           #227    ICSD_77473    binary_carbide_ICSD_frequency=1
1284.AB/BA     hP9      P-31m                    C1W2            #162    ICSD_162419   binary_carbide_ICSD_frequency=1
1285.AB/BA     oP12     Pbcn                     C1W2            #60     ICSD_43017    binary_carbide_ICSD_frequency=1
1286.AB/BA     mP12     P2_{1}/c                 C2Ca1           #14     ICSD_24074    binary_carbide_ICSD_frequency=1
1287.AB/BA     cP12     Pn-3m                    C2Ca1           #224    ICSD_31092    binary_carbide_ICSD_frequency=1
1288.AB/BA     mS12     C2/m                     C2Ca1           #12     ICSD_54185    binary_carbide_ICSD_frequency=1
1289.AB/BA     aP12     P-1                      C2Ca1           #2      ICSD_66663    binary_carbide_ICSD_frequency=1
1290.AB/BA     oP20     Pnma                     C2Cr3           #62     ICSD_15086    binary_carbide_ICSD_frequency=1
1291.AB/BA     oS20     Cmcm                     C2Cr3           #63     ICSD_57009    binary_carbide_ICSD_frequency=1
1292.AB/BA     oP20     Pnma                     C2Cr3           #62     ICSD_653958   binary_carbide_ICSD_frequency=1
1293.AB/BA     oP16     Pnma                     C2Cs2           #62     ICSD_51531    binary_carbide_ICSD_frequency=1
1294.AB/BA     hP12     P-62m                    C2Cs2           #189    ICSD_51533    binary_carbide_ICSD_frequency=1
1295.AB/BA     mS12     C2/c                     C2Eu1           #15     ICSD_380495   binary_carbide_ICSD_frequency=1
1296.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_423888   binary_carbide_ICSD_frequency=1
1297.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_89328    binary_carbide_ICSD_frequency=1
1298.AB/BA     mS12     C2/m                     C2Ir1           #12     ICSD_181488   binary_carbide_ICSD_frequency=1
1299.AB/BA     tI12     I4/mcm                   C2Ir1           #140    ICSD_181489   binary_carbide_ICSD_frequency=1
1300.AB/BA     tI6      I4/mmm                   C2Ir1           #139    ICSD_181490   binary_carbide_ICSD_frequency=1
1301.AB/BA     tI32     I4_{1}/acd               C2K2            #142    ICSD_36142    binary_carbide_ICSD_frequency=1
1302.AB/BA     tI64     I4_{1}/acd               C2K2            #142    ICSD_89528    binary_carbide_ICSD_frequency=1
1303.AB/BA     oI8      Immm                     C2Li2           #71     ICSD_25705    binary_carbide_ICSD_frequency=1
1304.AB/BA     tP6      P4_{2}/mnm               C2Mg1           #136    ICSD_88057    binary_carbide_ICSD_frequency=1
1305.AB/BA     mS28     C2/c                     C2Mn5           #15     ICSD_618246   binary_carbide_ICSD_frequency=1
1306.AB/BA     mS28     C2/c                     C2Mn5           #15     ICSD_69535    binary_carbide_ICSD_frequency=1
1307.AB/BA     cP12     P4_{3}32                 C2Na1           #212    ICSD_186177   binary_carbide_ICSD_frequency=1
1308.AB/BA     tI32     I4_{1}/acd               C2Na2           #142    ICSD_28066    binary_carbide_ICSD_frequency=1
1309.AB/BA     oI8      Immm                     C2Na2           #71     ICSD_95835    binary_carbide_ICSD_frequency=1
1310.AB/BA     hP6      P6_{3}/mmc               C2Re1           #194    ICSD_184660   binary_carbide_ICSD_frequency=1
1311.AB/BA     tI6      I4/mmm                   C2Re1           #139    ICSD_184661   binary_carbide_ICSD_frequency=1
1312.AB/BA     tP3      P4/mmm                   C2Re1           #123    ICSD_184664   binary_carbide_ICSD_frequency=1
1313.AB/BA     tP6      P4_{2}/mmc               C2Si1           #131    ICSD_187721   binary_carbide_ICSD_frequency=1
1314.AB/BA     mS12     C2/c                     C2Th1           #15     ICSD_26568    binary_carbide_ICSD_frequency=1
1315.AB/BA     tP6      P4_{2}/mmc               C2Th1           #131    ICSD_77469    binary_carbide_ICSD_frequency=1
1316.AB/BA     hP20     P6_{3}mc                 C3Fe7           #186    ICSD_76830    binary_carbide_ICSD_frequency=1
1317.AB/BA     oP16     Pnma                     C3Ir1           #62     ICSD_181492   binary_carbide_ICSD_frequency=1
1318.AB/BA     oP16     Pbam                     C3Ir5           #55     ICSD_181485   binary_carbide_ICSD_frequency=1
1319.AB/BA     oP10     Pnnm                     C3Mg2           #58     ICSD_71941    binary_carbide_ICSD_frequency=1
1320.AB/BA     cP7      Pm-3m                    C3Nb4           #221    ICSD_42758    binary_carbide_ICSD_frequency=1
1321.AB/BA     tP20     P-4c2                    C3Os2           #116    ICSD_168275   binary_carbide_ICSD_frequency=1
1322.AB/BA     oP40     Pbcn                     C3Os2           #60     ICSD_168276   binary_carbide_ICSD_frequency=1
1323.AB/BA     cI28     I-43d                    C3Sc4           #220    ICSD_42760    binary_carbide_ICSD_frequency=1
1324.AB/BA     mP10     P2_{1}/c                 C4Ir1           #14     ICSD_181493   binary_carbide_ICSD_frequency=1
1325.AB/BA     oF20     Ccce                     C4Ir1           #68     ICSD_181496   binary_carbide_ICSD_frequency=1
1326.AB/BA     tP5      P4/nbm                   C4Ir1           #125    ICSD_181497   binary_carbide_ICSD_frequency=1
1327.AB/BA     oS40     P6_{3}/mmc               C4Ir1           #194    ICSD_181498   binary_carbide_ICSD_frequency=1
1328.AB/BA     tP36     P-4c2                    C5Ir4           #116    ICSD_181486   binary_carbide_ICSD_frequency=1
1329.AB/BA     mS22     C2/m                     C5Nb6           #12     ICSD_20695    binary_carbide_ICSD_frequency=1
1330.AB/BA     hR13     R-3m                     C5Ti8           #166    ICSD_20822    binary_carbide_ICSD_frequency=1
1331.AB/BA     hP33     P3_{1}12                 C5V6            #151    ICSD_71098    binary_carbide_ICSD_frequency=1
1332.AB/BA     cP126    Pm-3                     C60Sr3          #200    ICSD_75354    binary_carbide_ICSD_frequency=1
1333.AB/BA     hP14     P6_{3}/mmc               C6Eu1           #194    ICSD_169041   binary_carbide_ICSD_frequency=1
1334.AB/BA     hP14     P6_{3}/mmc               C6Yb1           #194    ICSD_601565   binary_carbide_ICSD_frequency=1
1335.AB/BA     mP44     P2_{1}/c                 C7Ho4           #14     ICSD_74817    binary_carbide_ICSD_frequency=1
1336.AB/BA     cI40     Im-3m                    C7Ir3           #229    ICSD_181491   binary_carbide_ICSD_frequency=1
1337.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_22177    binary_carbide_ICSD_frequency=1
1338.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_43259    binary_carbide_ICSD_frequency=1
1339.AB/BA     hP9      P6/mmm                   C8Cs1           #191    ICSD_74641    binary_carbide_ICSD_frequency=1
1340.AB/BA     oF72     Fddd                     C8K1            #70     ICSD_70020    binary_carbide_ICSD_frequency=1
1341.AB/BA     oF72     Fddd                     C8Rb1           #70     ICSD_200563   binary_carbide_ICSD_frequency=1


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

ANRL Label                              Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
*****************************************************************************************************************************************************
ANRL - BCC-SQS
*****************************************************************************************************************************************************
AB_aP16_2_4i_4i.AB                      aP16      P-1                 -         AB                  #2        bcc-SQS [10.1103/PhysRevB.69.214202] (O. Levy)
A5B11_mP16_6_2abc_2a3b3c.AB/BA          mP16      Pm                  -         A5B11               #6        bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)      
AB3_mC32_8_4a_12a.AB/BA                 mC32      Cm                  -         AB3                 #8        bcc-SQS [10.1103/PhysRevB.69.214202] (O. Levy)
AB3_mC32_8_4a_4a4b.AB/BA                mC32      Cm                  -         AB3                 #8        bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)
A3B13_oC32_38_ac_a2bcdef.AB/BA          oC32      Amm2                -         A3B13               #38       bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)
A3B5_oC32_38_abce_abcdf.AB/BA           oC32      Amm2                -         A3B5                #38       bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)
AB7_hR16_166_c_c2h.AB/BA                hR16      R-3m                -         AB7                 #166      bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)

****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2021                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
