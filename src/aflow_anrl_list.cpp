// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Stefano Curtarolo, David Hicks - 2016

#ifndef _AFLOW_ANRL_LIST_CPP // AFLOW_REMOVE_GREP
#define _AFLOW_ANRL_LIST_CPP // AFLOW_REMOVE_GREP
#include "aflow.h" // AFLOW_REMOVE_GREP  //CO20200521

// ***************************************************************************
namespace anrl { 
  uint PrototypeANRL_LoadList(vector<string>& vproto,
      vector<string>& vproto_label,
      vector<uint>& vproto_nspecies,
      vector<uint>& vproto_natoms,
      vector<uint>& vproto_spacegroup,
      vector<uint>& vproto_nunderscores,
      vector<uint>& vproto_nparameters,
      vector<string>& vproto_Pearson_symbol,
      vector<string>& vproto_params,
      vector<string>& vproto_Strukturbericht,
      vector<string>& vproto_prototype,
      vector<string>& vproto_dialect) {

    vproto.clear();
    vproto_label.clear();
    vproto_nspecies.clear();
    vproto_natoms.clear();
    vproto_spacegroup.clear();
    vproto_nunderscores.clear();
    vproto_nparameters.clear();
    vproto_Pearson_symbol.clear();
    vproto_params.clear();
    vproto_Strukturbericht.clear();
    vproto_prototype.clear();
    vproto_dialect.clear();


    //Label     # of Species    # atoms/primcell    #space_group_number   # Underscores     # Parameters    pearson_symbol    params    Strukturbericht     Prototype   Dialect        
    // -------------------------------------------------------------------------
    // Part 1
    // -------------------------------------------------------------------------
    vproto.push_back("AB2_aP12_1_4a_8a;2;12;1;4;42;aP12;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;FeS2;FeS2");
    vproto.push_back("ABC2_aP16_1_4a_4a_8a;3;16;1;5;54;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;AsKSe2;AsKSe2");
    vproto.push_back("A2B_aP6_2_2i_i;2;6;2;4;15;aP6;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;P2I4;P2I4");
    vproto.push_back("A_aP4_2_aci;1;4;2;3;9;aP4;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3;-;Cf;Cf");
    vproto.push_back("A2B_mP12_3_bc3e_2e;2;12;3;4;21;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;SiO2;SiO2");
    vproto.push_back("A_mP4_4_2a;1;4;4;3;10;mP4;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2;-;Te;High-Pressure Te");
    vproto.push_back("A_mC12_5_3c;1;6;5;3;13;mC12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;A19;Po;Po");
    vproto.push_back("A3BC_mC10_8_ab_a_a;3;5;8;5;13;mC10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,y4,z4;-;Pb(Zr0.52Ti0.48)O3;Monoclinic PZT [PbO3]");
    vproto.push_back("A2B_mC144_9_24a_12a;2;72;9;4;112;mC144;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36;-;SiO2;Monoclinic Low Tridymite");
    vproto.push_back("AB_mP4_11_e_e;2;4;11;4;8;mP4;a,b/a,c/a,beta,x1,z1,x2,z2;-;NiTi;NiTi2");
    vproto.push_back("ABC3_mP10_11_e_e_ef;3;10;11;5;13;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,y4,z4;G0_6;KClO3;KClO3");
    vproto.push_back("A_mP16_11_8e;1;16;11;3;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;alphaPu;alphaPu");
    vproto.push_back("AB2_mC6_12_a_i;2;3;12;4;6;mC6;a,b/a,c/a,beta,x2,z2;C34/-/-/-/-;AuTe2/Hg1O2/Ni1O2;Calaverite/Hg1O2 (ICSD #48214)/Ni1O2 (ICSD #88720)"); //DX20210105 - added two metal-oxide prototypes
    vproto.push_back("A_mC34_12_ah3i2j;1;17;12;3;17;mC34;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;betaPu;betaPu");
    vproto.push_back("AB3_mC16_12_g_ij;2;8;12;4;10;mC16;a,b/a,c/a,beta,y1,x2,z2,x3,y3,z3;D0_15;AlCl3;AlCl3");
    vproto.push_back("A5B2_mC14_12_a2i_i;2;7;12;4;10;mC14;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4;-/-;Au5Mn2/B5W2;Au5Mn2/B5W2 (ICSD #20326)"); //DX20210104 - added metal-boride
    vproto.push_back("A_mC4_12_i;1;2;12;3;6;mC4;a,b/a,c/a,beta,x1,z1;-;O2;alphaO2");
    vproto.push_back("ABC4_mP12_13_e_a_2g;3;12;13;5;11;mP12;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;E1_b;AgAuTe4;Sylvanite");
    vproto.push_back("A_mP84_13_21g;1;84;13;3;67;mP84;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;P;Monoclinic Phosphorus");
    vproto.push_back("A2B_mP12_14_2e_e;2;12;14;4;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;C43/-/-/-/-/-/-/-/-/-;ZrO2/O2W/O2V/N2Os1/O2Zr1/O2Ti1/O2Ti1/O2Tc1/O2Zr1/C2Ca1;Baddeleyite/O2W (ICSD #8217)/O2V (ICSD #647610)/N2Os1 (ICSD #240759)/O2Zr1 (ICSD #15983)/O2Ti1 (ICSD #154035)/O2Ti1 (ICSD #154036)/O2Tc1 (ICSD #173153)/O2Zr1 (ICSD #659226)/C2Ca1 (ICSD #24074)"); //DX20200703 - added R. Friedrich binary oxide //DX20201019 - added metal nitride //DX20210105 - added five metal-oxide prototypes //DX20210120 - added metal-carbide
    vproto.push_back("A_mP32_14_8e;1;32;14;3;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;A_l;Se;betaSe");
    vproto.push_back("A_mP64_14_16e;1;64;14;3;52;mP64;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;A_k;Se;Se");
    vproto.push_back("A2B5_mC28_15_f_e2f;2;14;15;4;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-;B2Pd5/Sb2O5/Nb2O5;B2Pd5/Sb2O5 (ICSD #1422)/zeta-Nb2O5 (B-Nb2O5)"); //DX20200703 - added R. Friedrich binary oxide //DX20210427 - added zeta-Nb2O5 (B-Nb2O5) from part 3
    vproto.push_back("AB_mC8_15_c_e;2;4;15;4;5;mC8;a,b/a,c/a,beta,y2;B26/-;CuO/Ag1O1;Tenorite/Ag1O1 (ICSD #27659)"); //DX20210105 - added metal-oxide prototype
    vproto.push_back("A2B_mC48_15_ae3f_2f;2;24;15;4;20;mC48;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;SiO2;Coesite");
    vproto.push_back("ABC6D2_mC40_15_e_e_3f_f;4;20;15;6;18;mC40;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;- (S4_{1});CaFeO6Si2 (CaMgO6Si2);Esseneite and Diopside (CaMg(SiO3)2, S4_{1})"); //DX20210428 - added equivalent part 3 prototype (diopside, CaMg(SiO3)2, S4_{1}, http://aflow.org/prototype-encyclopedia/ABC6D2_mC40_15_e_e_3f_f.S4_1.html)
    vproto.push_back("ABC4_oP12_16_ag_cd_2u;3;12;16;5;9;oP12;a,b/a,c/a,x5,y5,z5,x6,y6,z6;-;AlPS4;AlPS4");
    vproto.push_back("AB3_oP16_18_ab_3c;2;16;18;4;14;oP16;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;D0_{17};BaS3;BaS3"); //DX 20180925 - moved Strukturberict designation from AB3_tP8_113_a_ce per M. Mehl's suggestion (historical reasons)
    vproto.push_back("A2B_oP12_19_2a_a;2;12;19;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ag2Se;Naumannite");
    vproto.push_back("A2B_oC24_20_abc_c;2;12;20;4;11;oC24;a,b/a,c/a,x1,y2,x3,y3,z3,x4,y4,z4;-;SiO2;Orthorhombic Tridymite");
    vproto.push_back("AB_oP2_25_b_a;2;2;25;4;5;oP2;a,b/a,c/a,z1,z2;-;CdTe;High-Pressure CdTe");
    vproto.push_back("AB2_oP24_28_acd_2c3d;2;24;28;4;22;oP24;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;C46;AuTe2;Krennerite");
    vproto.push_back("AB3C4_oP16_31_a_ab_2ab;3;16;31;5;17;oP16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;H2_5;AsCu3S4;Enargite");
    vproto.push_back("AB_oP8_33_a_a;2;8;33;4;9;oP8;a,b/a,c/a,x1,y1,z1,x2,y2,z2;-;CoAs;Modderite");
    vproto.push_back("AB3C4_oP32_33_a_3a_4a;3;32;33;5;27;oP32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;AsK3S4;AsK3S4");
    vproto.push_back("A2B_oC12_36_2a_a;2;6;36;4;9;oC12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;C24;HgBr2;HgBr2");
    vproto.push_back("A2BC_oC8_38_e_a_b;3;4;38;5;7;oC8;a,b/a,c/a,z1,z2,y3,z3;-;C2CeNi;C2CeNi");
    vproto.push_back("A2B_oC12_38_de_ab;2;6;38;4;9;oC12;a,b/a,c/a,z1,z2,y3,z3,y4,z4;-;Au2V;Au2V");
    vproto.push_back("AB4_oC20_41_a_2b;2;10;41;4;10;oC20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;D1_c;PtSn4;PtSn4");
    vproto.push_back("AB2_oC24_41_2a_2b;2;12;41;4;11;oC24;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4;C_e;PdSn2;PdSn2");
    vproto.push_back("AB2_oF72_43_ab_3b;2;18;43;4;16;oF72;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;C44;GeS2;GeS2");
    vproto.push_back("AB_oI4_44_a_b;2;2;44;4;5;oI4;a,b/a,c/a,z1,z2;-;AsGa;High-pressure GaAs");
    vproto.push_back("A2B3C7D_oP13_47_t_aq_eqrs_h;4;13;47;6;8;oP13;a,b/a,c/a,z4,z5,z6,z7,z8;-;YBa2Cu3O7-x;1212C [YBa2Cu3O7-x]");
    vproto.push_back("AB_oP4_51_e_f;2;4;51;4;5;oP4;a,b/a,c/a,z1,z2;B19;AuCd;beta'-AuCd");
    vproto.push_back("A3B2_oP20_56_ce_e;2;20;56;4;10;oP20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;D5_11;Sb2O3;Sb2O3");
    vproto.push_back("ABCD_oP16_57_d_c_d_d;4;16;57;6;10;oP16;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;F5_9;KCNS;KCNS");
    vproto.push_back("AB_oP8_57_d_d;2;8;57;4;7;oP8;a,b/a,c/a,x1,y1,x2,y2;-;TlF;TlF-II, 57");
    vproto.push_back("AB2_oP6_58_a_g;2;6;58;4;5;oP6;a,b/a,c/a,x2,y2;C35/-/C18/-/-;CaCl2/etaFe2C/FeS2/Ir1N2/Au1N2;Hydrophilite/eta-Fe2C/Marcasite/Ir1N2 (ICSD #160620)/Au1N2 (ICSD #166465)"); //DX20201019 - added two metal nitrides
    vproto.push_back("AB_oP4_59_a_b;2;4;59;4;5;oP4;a,b/a,c/a,z1,z2;-;CuTe;Vulcanite");
    vproto.push_back("ABC_oP6_59_a_a_a;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;-;CNCl;CNCl");
    vproto.push_back("A3B_oP8_59_bf_a;2;8;59;4;7;oP8;a,b/a,c/a,z1,z2,x3,z3;D0_a;TiCu3;betaTiCu3");
    vproto.push_back("AB_oP16_61_c_c;2;16;61;4;9;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2;B_e/-;CdSb/C2N2;CdSb/C2N2 (ICSD #15870)"); //DX20210219 - added carbo-nitride
    vproto.push_back("A2B_oP24_61_2c_c;2;24;61;4;12;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;C21/C52;TiO2/O2Te;Brookite/Tellurite (beta-TeO2, C52)"); //DX20210427 - added Tellurite (beta-TeO2, C52) from part 3
    vproto.push_back("A3B2_oP20_62_3c_2c;2;20;62;4;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;D5_8;Sb2S3;Stibnite");
    vproto.push_back("AB3C_oP20_62_c_cd_a;3;20;62;5;10;oP20;a,b/a,c/a,x2,z2,x3,z3,x4,y4,z4;-/-;CaTiO3/CaO3Zr;CaTiO3 Pnma Perovskite/CaO3Zr (ICSD #97463)"); //DX20200703 - added R. Friedrich ternary oxide
    vproto.push_back("A4B_oP20_62_2cd_c;2;20;62;4;12;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4;-;MgB4;MgB4");
    vproto.push_back("AB2C_oP16_62_c_2c_c;3;16;62;5;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;F5_6/-/-;CuSbS2/C1N2Sr1/C1N2Pb1;Chalcostibite/C1N2Sr1 (ICSD #75040)/C1N2Pb1 (ICSD #410915)"); //DX20210221 - added two carbo-nitrides
    vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C37/C25/C23/C29/-/-/-/-/-/C53;Co2Si/HgCl2/PbCl2/SrH2/N2Re1/O2Zr1/O2Sn1/O2Ti1/O2Pr1/Br2Sr;Co2Si/HgCl2/Cotunnite/SrH2/N2Re1 (ICSD #187446)/O2Zr1 (ICSD #56696)/O2Sn1 (ICSD #181282)/O2Ti1 (ICSD #182577)/O2Pr1 (ICSD #380398)/C53 (SrBr2) (obsolete)"); //added info from part 2 A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C29;SrH2;SrH2 //DX20201019 - added metal-nitride //DX20210106 - added four metal-oxide //DX20210427 - added C53 (SrBr2) (obsolete) from part 3
    vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B16/B31/B27/B29/B14 (B_{d})/-/-/-/-;GeS/MnP/FeB/SnS/FeAs (NiSi)/oxide/B1Ti1/B1Li1/B1Co1;GeS/MnP/FeB/SnS/Westerveldite and eta-NiSi/oxide/B1Ti1 (ICSD #24701)/B1Li1;B1Li1 (ICSD #153291)/B1Co1 (ICSD #612863)"); //added info from part 2 AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B14;FeAs;Westerveldite //DX20210104 - added three metal-borides //DX20210428 - added equivalent part 3 prototype (eta-NiSi, http://aflow.org/prototype-encyclopedia/AB_oP8_62_c_c.NiSi.html)
    vproto.push_back("AB3_oP16_62_c_cd;2;16;62;4;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;D0_11/-/-/-;Fe3C/N1Na3/C1Fe3/C1Cr3;Cementite/N1Na3 (ICSD #165990)/C1Fe3 (ICSD #167667)/C1Cr3 (ICSD #617486)"); //DX20201019 - added metal-nitride //DX20210120 - added two metal-carbides
    vproto.push_back("A3B7_oP40_62_cd_3c2d;2;40;62;4;20;oP40;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;D10_1;C3Cr7;C3Cr7");
    vproto.push_back("A_oP8_62_2c;1;8;62;3;7;oP8;a,b/a,c/a,x1,z1,x2,z2;A_c;alphaNp;alphaNp");
    vproto.push_back("AB2C_oC16_63_c_2c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-/-/-;SrCuO2/Au1Ca2N1/Co1Sn2Tb1;SrCuO2/Au1Ca2N1 (ICSD #85528)/Co1Sn2Tb1 (ICSD #240096)"); //DX20201016 - added metal-nitride //DX20201028 - added one metal
    vproto.push_back("A2B_oC12_63_2c_c;2;6;63;4;6;oC12;a,b/a,c/a,y1,y2,y3;C49;ZrSi2;ZrSi2");
    vproto.push_back("AB_oC8_63_c_c;2;4;63;4;5;oC8;a,b/a,c/a,y1,y2;B33/-/-;CrB/-/Al1N1;CrB/-/Al1N1 (ICSD #163951)"); //DX20201019 - added metal-nitride
    vproto.push_back("A_oC4_63_c;1;2;63;3;4;oC4;a,b/a,c/a,y1;A20;alpha-U;alpha-Uranium");
    vproto.push_back("A_oC8_64_f;1;4;64;3;5;oC8;a,b/a,c/a,y1,z1;A11/A17/A14;alphaGa/P/I2;alphaGallium/Black Phosphorus/Molecular Iodine");
    vproto.push_back("A2B2C_oC80_64_efg_efg_df;3;40;64;5;18;oC80;a,b/a,c/a,x1,y2,y3,y4,z4,y5,z5,y6,z6,x7,y7,z7,x8,y8,z8;-;MgB2C2;MgB2C2");
    vproto.push_back("AB_oC8_65_j_g;2;4;65;4;5;oC8;a,b/a,c/a,x1,y2;-;alphaIrV;alphaIrV");
    vproto.push_back("A3B5_oC16_65_ah_bej;2;8;65;4;5;oC16;a,b/a,c/a,x4,y5;-;Ga3Pt5;Ga3Pt5");
    vproto.push_back("AB3_oC8_65_a_bf;2;4;65;4;3;oC8;a,b/a,c/a;L1_3;CdPt3;Predicted CdPt3");
    vproto.push_back("AB_oF8_69_a_b;2;2;69;4;3;oF8;a,b/a,c/a;B24;TlF;TlF");
    vproto.push_back("A_oF8_70_a;1;2;70;3;3;oF8;a,b/a,c/a;-;gammaPu;gamma-Pu");
    vproto.push_back("A2B_oF24_70_e_a;2;6;70;4;4;oF24;a,b/a,c/a,x2;C54;TiSi2;TiSi2");
    vproto.push_back("A_oF128_70_4h;1;32;70;3;15;oF128;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;A16;S;alpha-Sulfur");
    vproto.push_back("AB2_oI6_71_a_i;2;3;71;4;4;oI6;a,b/a,c/a,z2;-;ReSi2;ReSi2");
    vproto.push_back("AB2_oI6_71_a_g;2;3;71;4;4;oI6;a,b/a,c/a,y2;-;MoPt2;MoPt2");
    vproto.push_back("A2B_oI12_72_j_a;2;6;72;4;5;oI12;a,b/a,c/a,x2,y2;C42;SiS2;SiS2");
    vproto.push_back("AB4C_tI12_82_c_g_a;3;6;82;5;5;tI12;a,c/a,x3,y3,z3;H0_7;BPO4;BPO4");
    vproto.push_back("A2BC4_tI14_82_bc_a_g;3;7;82;5;5;tI14;a,c/a,x4,y4,z4;E3;CdAl2S4;CdAl2S4");
    vproto.push_back("AB_tP16_84_cej_k;2;16;84;4;7;tP16;a,c/a,x3,y3,x4,y4,z4;B34;PdS;PdS");
    vproto.push_back("A4B5_tI18_87_h_ah;2;9;87;4;6;tI18;a,c/a,x2,y2,x3,y3;-;Ti5Te4;Ti5Te4");
    vproto.push_back("AB4_tI10_87_a_h;2;5;87;4;4;tI10;a,c/a,x2,y2;D1_a/-;Ni4Mo/C1Fe4;Ni4Mo/C1Fe4 (ICSD #187144)"); //DX20210120 - added metal-carbide
    vproto.push_back("A2B_tP12_92_b_a;2;12;92;4;6;tP12;a,c/a,x1,x2,y2,z2;-/-;SiO2/O2Ti (O2Te);alpha Cristobalite/O2Ti and alpha-O2Te (paratellurite)"); //DX20210428 - added equivalent part 3 prototype (alpha-TeO2, paratellurite, http://aflow.org/prototype-encyclopedia/A2B_tP12_92_b_a.TeO2.html)
    vproto.push_back("A2B_tP36_96_3b_ab;2;36;96;4;15;tP36;a,c/a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;SiO2;Keatite");
    vproto.push_back("A_tP12_96_ab;1;12;96;3;6;tP12;a,c/a,x1,x2,y2,z2;-;Si;\"ST12'' of Si");
    vproto.push_back("A3BC_tP5_99_bc_a_b;3;5;99;5;6;tP5;a,c/a,z1,z2,z3,z4;-;Pb(Zr0.52Ti0.48)O3;Tetragonal PZT [PbO3]");
    vproto.push_back("AB3_tP8_113_a_ce;2;8;113;4;5;tP8;a,c/a,z2,x3,z3;-;BaS3;BaS3");//DX 20180925 - moved Strukturberict designation to AB3_oP16_18_ab_3c per M. Mehl's suggestion (historical reasons)
    vproto.push_back("A2BC4D_tI16_121_d_a_i_b;4;8;121;6;4;tI16;a,c/a,x4,z4;H2_6;Cu2FeS4Sn;Stannite");
    vproto.push_back("ABC2_tI16_122_a_b_d;3;8;122;5;3;tI16;a,c/a,x3;E1_1/-;CuFeS2/C1Mg1N2;Chalcopyrite/C1Mg1N2 (ICSD #44110)"); //DX20210221 - added carbo-nitride
    vproto.push_back("AB5C_tP7_123_b_ci_a;3;7;123;5;3;tP7;a,c/a,z4;-;HoCoGa5;HoCoGa5");
    vproto.push_back("AB3_tP4_123_a_ce;2;4;123;4;2;tP4;a,c/a;L6_0;CuTi3;CuTi3");
    vproto.push_back("AB_tP2_123_a_d;2;2;123;4;2;tP2;a,c/a;L1_0/L2_{a};CuAu/CuTi;CuAu/delta-CuTi (L2a)"); //DX20210427 - added delta-CuTi (L2a) from part 3
    vproto.push_back("ABC2_tP4_123_d_a_f;3;4;123;5;2;tP4;a,c/a;-;CaCuO2;CaCuO2");
    vproto.push_back("A2B3_tP10_127_g_ah;2;10;127;4;4;tP10;a,c/a,x2,x3;D5_a;Si2U3;Si2U3");
    vproto.push_back("ABCD_tP8_129_c_b_a_c;4;8;129;6;4;tP8;a,c/a,z3,z4;-;AsCuSiZr;AsCuSiZr");
    vproto.push_back("A_tP4_129_ac;1;4;129;3;3;tP4;a,c/a,z2;A_d;betaNp;betaNp");
    vproto.push_back("ABC_tP6_129_c_a_c;3;6;129;5;4;tP6;a,c/a,z2,z3;E0_1/-;PbFCl/K1Na1O1;Matlockite/K1Na1O1 (ICSD #32743)"); //DX20210106 - added metal-oxide
    vproto.push_back("A2B_tP6_129_ac_c;2;6;129;4;4;tP6;a,c/a,z2,z3;C38;Cu2Sb;Cu2Sb");
    vproto.push_back("AB_tP4_129_a_c;2;4;129;4;3;tP4;a,c/a,z2;B10;PbO;PbO");
    vproto.push_back("AB_tP4_129_c_c;2;4;129;4;4;tP4;a,c/a,z1,z2;B11/-;gammaCuTi/-;gammaCuTi/-");
    vproto.push_back("AB_tP4_131_c_e;2;4;131;4;2;tP4;a,c/a;B17/-;PtS/PdO;PtS/PdO (ICSD #26598)"); //DX20200703 - added R. Friedrich binary oxide
    vproto.push_back("A_tP50_134_b2m2n;1;50;134;3;12;tP50;a,c/a,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5;A_g;B;T-50 Boron");
    vproto.push_back("A_tP30_136_bf2ij;1;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;A_b;betaU;betaU");
    vproto.push_back("AB_tP8_136_g_f;2;8;136;4;4;tP8;a,c/a,x1,x2;-;betaBeO;betaBeO");
    vproto.push_back("A2B_tP6_136_f_a;2;6;136;4;3;tP6;a,c/a,x2;C4/-/-;TiO2/O2W1/C2Mg1;CrSi2/O2W1 (ICSD #647647)/C2Mg1 (ICSD #88057)"); //DX20210106 - added metal-oxide //DX20210120 - added metal-carbide
    //vproto.push_back("sigma_tP30_136_bf2ij;5;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;D8_b;CrFe;sigmaCrFe");
    vproto.push_back("sigma_tP30_136_bf2ij;1;30;136;3;9;tP30;a,c/a,x2,x3,y3,x4,y4,x5,z5;D8_b;CrFe;sigmaCrFe");
    vproto.push_back("A_tP4_136_f;1;4;136;3;3;tP4;a,c/a,x1;-;N2;gammaN2");
    vproto.push_back("A_tP16_138_j;1;16;138;3;5;tP16;a,c/a,x1,y1,z1;A18;Cl2;Cl2");
    vproto.push_back("A3B_tI16_139_cde_e;2;8;139;4;4;tI16;a,c/a,z3,z4;D0_23;Al3Zr;Al3Zr");
    vproto.push_back("A_tI4_139_e;1;2;139;3;3;tI4;a,c/a,z1;-;Si;Hypothetical BCT5 Si");
    vproto.push_back("AB2C4_tI14_139_a_e_ce;3;7;139;5;4;tI14;a,c/a,z3,z4;-/-/-;(La,Ba)2CuO4/O4Sr2Ti (F4K2Ni)/Cu1Nd2O4;0201 [(La,Ba)2CuO4]/O4Sr2Ti (ICSD #157402) and K2NiF4 (part 3)/Cu1Nd2O4 (ICSD #86753)"); //DX20200703 - added R. Friedrich ternary oxide //DX20210428 - added equivalent part 3 prototype (K2NiF4, http://aflow.org/prototype-encyclopedia/A4B2C_tI14_139_ce_e_a.html)
    vproto.push_back("A12B_tI26_139_fij_a;2;13;139;4;4;tI26;a,c/a,x3,x4;D2_b;Mn12Th;Mn12Th");
    vproto.push_back("A_tI2_139_a;1;1;139;3;2;tI2;a,c/a;A6/A_a;In/Pa;In/alphaPa");
    vproto.push_back("A_tI8_139_h;1;4;139;3;3;tI8;a,c/a,x1;-;C;Hypothetical Tetrahedrally Bonded Carbon with 4-Member Rings");
    vproto.push_back("A3B_tI8_139_bd_a;2;4;139;4;2;tI8;a,c/a;D0_22;Al3Ti;Al3Ti");
    vproto.push_back("AB2_tI6_139_a_e;2;3;139;4;3;tI6;a,c/a,z2;C11_b/-/-/C11_{a}/-/-;MoSi2/-/-/KO2 (C2Ca)/Ca1O2/Ba1O2;MoSi2/-/(ICSD #24248)/KO2 (ICSD #38245) and CaC2-I (part 3)/Ca1O2 (ICSD #20275)/Ba1O2 (ICSD #24729)"); //DX20200703 - added R. Friedrich binary oxide //DX20210106 - added two metal-oxide prototypes //DX20210428 - added equivalent part 3 prototype (C2Ca-I, C11_{a}, http://aflow.org/prototype-encyclopedia/A2B_tI6_139_e_a.html, part 3)
    vproto.push_back("A4B5_tI18_139_i_ah;2;9;139;4;4;tI18;a,c/a,x2,x3;-;V4Zn5;V4Zn5");
    vproto.push_back("A4B_tI10_139_de_a;2;5;139;4;3;tI10;a,c/a,z3;D1_3;Al4Ba;Al4Ba");
    vproto.push_back("A8B_tI18_139_hi_a;2;9;139;4;4;tI18;a,c/a,x2,x3;-;Pt8Ti;Pt8Ti");
    vproto.push_back("A2B_tI6_139_d_a;2;3;139;4;2;tI6;a,c/a;L\'2;ThH2;ThH2");
    vproto.push_back("A2B_tI12_140_h_a;2;6;140;4;3;tI12;a,c/a,x2;C16/-;Al2Cu/C2Ir1;Khatyrkite/C2Ir1 (ICSD #181489)"); //DX20210120 - added metal-carbide
    vproto.push_back("AB3_tI16_140_b_ah;2;8;140;4;3;tI16;a,c/a,x3;D0_c;SiU3;SiU3");
    vproto.push_back("AB_tI16_140_ab_h;2;8;140;4;3;tI16;a,c/a,x3;B37;SeTl;SeTl");
    vproto.push_back("A4BC_tI24_141_h_b_a;3;12;141;5;4;tI24;a,c/a,y3,z3;-;ZrSiO4;Zircon");
    vproto.push_back("A_tI4_141_a;1;2;141;3;2;tI4;a,c/a;A5;Sn;betaSn");
    vproto.push_back("A3B4_tI28_141_ad_h;2;14;141;4;4;tI28;a,c/a,y3,z3;-/-;Mn3O4/Mn3O4;Hausmannite/Mn3O4 (ICSD #643198)"); //DX20210106 - added metal-oxide
    vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C5/C_{c}/-/-;TiO2/alpha-ThSi2/Mo2N1/O2Ti1;Anatase/alpha-ThSi2/Mo2N1 (ICSD #30593)/O2Ti1 (ICSD #200392)"); // added info from part 2 - A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C_{c};alpha-ThSi2;alpha-ThSi2 //DX20201103 - added metal-nitride //DX20210106 - added metal-oxide
    vproto.push_back("AB_tI16_141_e_e;2;8;141;4;4;tI16;a,c/a,z1,z2;B_g;MoB;MoB");
    vproto.push_back("A2B_tI24_141_2e_e;2;12;141;4;5;tI28;a,c/a,z1,z2,z3;-;Ga2Hf;Ga2Hf");
    vproto.push_back("AB_tI8_141_a_b;2;4;141;4;2;tI8;a,c/a;\"40\"/-;NbP/-;NbP/-");
    vproto.push_back("A2B3_tI80_141_ceh_3h;2;40;141;4;11;tI80;a,c/a,z2,y3,z3,y4,z4,y5,z5,y6,z6;-;In2S3;betaIn2S3");
    vproto.push_back("ABC4_tI96_142_e_ab_2g;3;48;142;5;9;tI96;a,c/a,x3,x4,y4,z4,x5,y5,z5;-;PPrS4;PPrS4");
    vproto.push_back("A2B_hP9_147_g_ad;2;9;147;4;6;hP9;a,c/a,z2,x3,y3,z3;B_b;AgZn;zetaAgZn");
    vproto.push_back("AB_hR16_148_cf_cf;2;16;148;4;10;hR16;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4;-;C8H8;Solid Cubane");
    vproto.push_back("AB3_hR8_148_c_f;2;8;148;4;6;hR8;a,c/a,x1,x2,y2,z2;D0_5;BiI3;BiI3");
    vproto.push_back("AB_hR26_148_b2f_a2f;2;26;148;4;14;hR26;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;PdAl;PdAl");
    vproto.push_back("AB3C_hR10_148_c_f_c;3;10;148;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-/-/-;FeTiO3/Cd1O3Ti1/Fe1O3Ti1;Ilmenite/Cd1O3Ti1 (ICSD #15989)/Fe1O3Ti1 (ICSD #153674)"); //DX20210106 - added two metal-oxides
    vproto.push_back("A2B_hP9_150_ef_bd;2;9;150;4;5;hP9;a,c/a,z2,x3,x4;C22;Fe2P;Original Fe2P");
    vproto.push_back("A3B_hP24_151_3c_2a;2;24;151;4;13;hP24;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;D0_4;CrCl3;CrCl3");
    vproto.push_back("A2B_hP9_152_c_a;2;9;152;4;6;hP9;a,c/a,x1,x2,y2,z2;-/-;SiO2/O2Ti1;alphaQuartz/O2Ti1 (ICSD #41493)"); //DX20210106 - added metal-oxide
    vproto.push_back("A_hP3_152_a;1;3;152;3;3;hP3;a,c/a,x1;A8;Se;gammaSe");
    vproto.push_back("AB_hP6_154_a_b;2;6;154;4;4;hP6;a,c/a,x1,x2;B9;HgS;Cinnabar");
    vproto.push_back("AB3_hR8_155_c_de;2;8;155;4;5;hR8;a,c/a,x1,y2,y3;D0_14;AlF3;AlF3");
    vproto.push_back("A3B2_hR5_155_e_c;2;5;155;4;4;hR5;a,c/a,x1,y2;D5_e;Ni3S2;Hazelwoodite");
    vproto.push_back("AB_hR6_160_b_b;2;6;160;4;6;hR6;a,c/a,x1,z1,x2,z2;B13;NiS;Millerite");
    vproto.push_back("AB_hR6_160_3a_3a;2;6;160;4;8;hR6;a,c/a,x1,x2,x3,x4,x5,x6;-;CSi;Moissanite 9R structure");
    vproto.push_back("ABC3_hR10_161_a_a_b;3;10;161;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-/-/-/-;LiNbO3/Na1Nb1O3/Li1Nb1O3/Ga1La1O3;Ferroelectric LiNbO3/Na1Nb1O3 (ICSD #9645)/Li1Nb1O3 (ICSD #28296)/Ga1La1O3 (ICSD #51036)"); //DX20210106 - added three metal-oxides
    vproto.push_back("AB2_hP9_162_ad_k;2;9;162;4;4;hP9;a,c/a,x3,z3;-;NV2;betaV2N");
    vproto.push_back("AB2CD2_hP36_163_h_i_bf_i;4;36;163;6;10;hP36;a,c/a,z2,x3,x4,y4,z4,x5,y5,z5;F5_10;KAg(CN)2;KAg2");
    vproto.push_back("A3B2_hP5_164_ad_d;2;5;164;4;4;hP5;a,c/a,z2,z3;D5_13/-/-/-/-;Al3Ni2/Ca3N2/Ca3N2/Ca3N2/O3Sc2;Al3Ni2/Ca3N2 (ICSD #162795)/Ca3N2 (ICSD #162796)/Ca3N2 (ICSD #169727)/O3Sc2 (ICSD #160203)"); //DX20201019 - added three metal-nitrides //DX20210106 - added metal-oxide
    vproto.push_back("AB2_hP3_164_a_d;2;3;164;4;3;hP3;a,c/a,z2;C6/-;CdI2/C1Sc2;omega Phase/C1Sc2 (ICSD #280743)"); //DX20210120 - added metal-carbide
    vproto.push_back("A3B_hP24_165_adg_f;2;24;165;4;7;hP24;a,c/a,z2,x3,x4,y4,z4;-;H3Ho;H3Ho");
    vproto.push_back("AB_hR2_166_a_b;2;2;166;4;2;hR2;a,c/a;L1_1/-/B22;CuPt/Fe1O1/KS;CuPt/Fe1O1 (ICSD #82236)/K(SH) (B22)"); //DX20210106 - added metal-oxide //DX20210427 - added K(SH) (B22) from part 3
    vproto.push_back("A_hR2_166_c;1;2;166;3;3;hR2;a,c/a,x1;A7/-/-;As/C/O2;alphaAs/Rhombohedral Graphite/betaO2");
    vproto.push_back("A_hR1_166_a;1;1;166;3;2;hR1;a,c/a;A_i/A10;Po/Hg;betaPo/alphaHg");
    vproto.push_back("A7B6_hR13_166_ah_3c;2;13;166;4;7;hR13;a,c/a,x2,x3,x4,x5,z5;D8_5;Fe7W6;Fe7W6 mu-phase");
    vproto.push_back("A_hR3_166_ac;1;3;166;3;3;hR3;a,c/a,x2;C19;Sm;alphaSm");
    vproto.push_back("A2B3_hR5_166_c_ac;2;5;166;4;4;hR5;a,c/a,x2,x3;C33;Bi2Te3;Bi2Te3");
    vproto.push_back("A5B2_hR7_166_a2c_c;2;7;166;4;5;hR7;a,c/a,x2,x3,x4;D8_i/-;B5Mo2/B5W2;Mo2B5/B5W2 (ICSD #20326)"); //DX20210104 - added metal-boride
    vproto.push_back("A_hR12_166_2h;1;12;166;3;6;hR12;a,c/a,x1,z1,x2,z2;-;B;alphaBoron");
    vproto.push_back("ABC2_hR4_166_a_b_c;3;4;166;5;3;hR4;a,c/a,x3;F5_1/-/-/-/-/-/-/-;CrNiS2/-/AlLiO2/CuFeO2 (CuFeO2)/Ag1Cr1O2/Cu1La1O2/Cu1Eu1O2/Cu1In1O2;Caswellsilverite/-/AlLiO2 (ICSD #28288)/CuFeO2 (ICSD #246912) and Rhombohedral Delafossite (CuFeO2, part 3)/Ag1Cr1O2 (ICSD #4149)/Cu1La1O2 (ICSD #18102)/Cu1Eu1O2 (ICSD #18106)/Cu1In1O2 (ICSD #55687)"); //DX20200703 - added R. Friedrich ternary oxides (2) //DX20210106 - added four metal-oxides //DX20210428 - added equivalent part 3 prototype (rhombohedral delafossite, CuFeO2, http://aflow.org/prototype-encyclopedia/ABC2_hR4_166_a_b_c.CuFeO2.html)
    vproto.push_back("A_hR105_166_bc9h4i;1;105;166;3;33;hR105;a,c/a,x2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;B;betaBoron");
    vproto.push_back("A6B_hR7_166_g_a;2;7;166;4;3;hR7;a,c/a,x2;-;CaC6;CaC6");
    vproto.push_back("ABC3_hR10_167_a_b_e;3;10;167;5;3;hR10;a,c/a,x3;-/G0_1;LiNbO3/CaCO3;Paraelectric LiNbO3/Calcite");
    vproto.push_back("A2B3_hR10_167_c_e;2;10;167;4;4;hR10;a,c/a,x1,x2;D5_1/-/-/-;Al2O3/O3V2/Al2O3/Al2O3;Corundum/O3V2 (ICSD #94768)/Al2O3 (ICSD #43732)/Al2O3 (ICSD #99783)"); //DX20200703 - added R. Friedrich binary oxide //DX20210106 - added two metal-oxides
    vproto.push_back("A2B_hP18_180_fi_bd;2;18;180;4;4;hP18;a,c/a,z3,x4;C_a;Mg2Ni;Mg2Ni");
    vproto.push_back("AB2_hP9_180_d_j;2;9;180;4;3;hP9;a,c/a,x2;C40;CrSi2;CrSi2");
    vproto.push_back("A2B_hP9_180_j_c;2;9;180;4;3;hP9;a,c/a,x2;C8;SiO2;beta-Quartz");
    vproto.push_back("AB3_hP8_182_c_g;2;8;182;4;3;hP8;a,c/a,x2;-/-;Fe3C/B1Fe3;Bainite/B1Fe3 (ICSD #184958)"); //DX20210104 - added metal-boride
    vproto.push_back("A_hP4_186_ab;1;4;186;3;4;hP4;a,c/a,z1,z2;-;C;Buckled Graphite");
    vproto.push_back("AB_hP8_186_ab_ab;2;8;186;4;6;hP8;a,c/a,z1,z2,z3,z4;B5;SiC;Moissanite-4H SiC");
    vproto.push_back("AB_hP4_186_b_b;2;4;186;4;4;hP4;a,c/a,z1,z2;B4/-/-/-/-;ZnS/Ga1N1/La1N1/Mo1N1/C1Sc1;Wurtzite/Ga1N1 (ICSD #153888)/La1N1 (ICSD #162195)/Mo1N1 (ICSD #168369)/C1Sc1 (ICSD #189087)"); //DX20201019 - added three metal-nitride //DX20210120 - added metal-carbide
    vproto.push_back("AB_hP12_186_a2b_a2b;2;12;186;4;8;hP12;a,c/a,z1,z2,z3,z4,z5,z6;B6;SiC;Moissanite-6H SiC");
    vproto.push_back("A5B3C_hP18_186_2a3b_2ab_b;3;18;186;5;11;hP18;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;E9_4;Al5C3N;Al5C3N");
    vproto.push_back("AB_hP4_186_b_a;2;4;186;4;4;hP4;a,c/a,z1,z2;B12/-;BN/Ga1N1;Original BN/Ga1N1 (ICSD #159250)"); //DX20201019 - added metal-nitride
    vproto.push_back("ABC_hP3_187_a_d_f;3;3;187;5;2;hP3;a,c/a;-;BaPtSb;BaPtSb");
    vproto.push_back("AB_hP2_187_d_a;2;2;187;4;2;hP2;a,c/a;B_h/-/-;WC/Ir1N1/Mg1O1;Tungsten Carbide/Ir1N1 (ICSD #186245)/Mg1O1 (ICSD #166273)"); //DX20201103 - added metal-nitride //DX20210106 - added metal-oxide
    vproto.push_back("A2B_hP9_189_fg_bc;2;9;189;4;4;hP9;a,c/a,x3,x4;C22;Fe2P;Revised Fe2P");
    vproto.push_back("AB4C_hP6_191_a_h_b;3;6;191;5;3;hP6;a,c/a,z3;-;AlB4Mg;AlB4Mg");
    vproto.push_back("AB5_hP6_191_a_cg;2;6;191;4;2;hP6;a,c/a;D2_d;CaCu5;CaCu5");
    vproto.push_back("A_hP1_191_a;1;1;191;3;2;hP1;a,c/a;A_f;gammaHgSn6-10;Simple Hexagonal Lattice");
    vproto.push_back("A3B_hP4_191_bc_a;2;4;191;4;2;hP4;a,c/a;-;Li3N;Li3 Ni");
    vproto.push_back("AB2_hP3_191_a_d;2;3;191;4;2;hP3;a,c/a;C32;AlB2;Hexagonal omega");
    vproto.push_back("A2B_hP6_191_h_e;2;6;191;4;4;hP6;a,c/a,z1,z2;C_h;Cu2Te;Cu2Te");
    vproto.push_back("AB_hP6_191_f_ad;2;6;191;4;2;hP6;a,c/a;B35/-;CoSn/N1Ta1;CoSn/N1Ta1 (ICSD #25659)"); //DX20201019 - added metal-nitride
    vproto.push_back("AB_hP8_194_ad_f;2;8;194;4;3;hP8;a,c/a,z3;B_i;AsTi;AsTi");
    vproto.push_back("A_hP6_194_h;1;6;194;3;3;hP6;a,c/a,x1;-;C;Hypothetical Tetrahedrally Bonded Carbon with 3-Member Rings");
    vproto.push_back("AB_hP12_194_af_bf;2;12;194;4;4;hP12;a,c/a,z3,z4;-;CMo;CMo structure");
    vproto.push_back("A_hP4_194_ac;1;4;194;3;2;hP4;a,c/a;A3\';La;alphaLa");
    vproto.push_back("AB3_hP8_194_c_bf;2;8;194;4;3;hP8;a,c/a,z3;D0_18;Na3As;Na3As");
    vproto.push_back("AB2_hP6_194_b_f;2;6;194;4;3;hP6;a,c/a,z2;-;CaIn2;CaIn2");
    vproto.push_back("AB_hP4_194_c_d;2;4;194;4;2;hP4;a,c/a;B_k/-;BN/Al1N1;BN/Al1N1 (ICSD #163950)"); //DX20201019 - added metal-nitride
    vproto.push_back("ABC2_hP8_194_d_a_f;3;8;194;5;3;hP8;a,c/a,z3;-;AlCCr2;AlCCr2");
    vproto.push_back("A3B_hP8_194_h_c;2;8;194;4;3;hP8;a,c/a,x2;D0_19;Ni3Sn;Ni3Sn");
    vproto.push_back("A_hP4_194_bc;1;4;194;3;2;hP4;a,c/a;A9;C;Hexagonal Graphite structure");
    vproto.push_back("AB2_hP6_194_c_f;2;6;194;4;3;hP6;a,c/a,z2;C7/-;MoS2/N1Re2;Molybdenite/N1Re2 (ICSD #169882)"); //DX20201103 - added metal-nitride
    vproto.push_back("A5B2_hP14_194_abdf_f;2;14;194;4;4;hP14;a,c/a,z4,z5;D8_h;W2B5;W2B5");
    vproto.push_back("AB2_hP12_194_f_ah;2;12;194;4;4;hP12;a,c/a,z2,x3;C14;MgZn2;MgZn2 Hexagonal Laves");
    vproto.push_back("ABC_hP6_194_c_d_a;3;6;194;5;2;hP6;a,c/a;-/-;LiBC/Al1Au1Ti1;LiBC/Al1Au1Ti1 (ICSD #57506)"); //DX20201028 - added metal ternary
    vproto.push_back("A_hP4_194_f;1;4;194;3;3;hP4;a,c/a,z1;-;C;Lonsdaleite");
    vproto.push_back("AB2_hP6_194_c_ad;2;6;194;4;2;hP6;a,c/a;B8_2;Ni2In;Ni2In");
    vproto.push_back("AB3C4_hP16_194_c_af_ef;3;16;194;5;5;hP16;a,c/a,z3,z4,z5;-;AlN3Ti4;AlN3Ti4");
    vproto.push_back("A_hP2_194_c;1;2;194;3;2;hP2;a,c/a;A3;Mg;Hexagonal Close Packed");
    vproto.push_back("AB2_hP24_194_ef_fgh;2;24;194;4;6;hP24;a,c/a,z1,z2,z3,x5;C36;MgNi2;MgNi2 Hexagonal Laves");
    vproto.push_back("AB_hP12_194_df_ce;2;12;194;4;4;hP12;a,c/a,z3,z4;B18;CuS;Covellite");
    vproto.push_back("AB_hP4_194_c_a;2;4;194;4;2;hP4;a,c/a;B8_1/-/- (L'3_{0})/-/-;NiAs/N1Nb1/N1Ta1 (Fe2N)/N1Re1/N1Nb1;NiAs/N1Nb1 (ICSD #76384)/N1Ta1 (ICSD #105123) and Fe2N (L'3_{0}, part 3)/N1Re1 (ICSD #181300)/N1Nb1 (ICSD #185559)"); //DX20201103 - added 4 metal-nitrides //DX20210428 - added equivalent part 3 prototype (Fe2N, L'3_{0}, http://aflow.org/prototype-encyclopedia/AB_hP4_194_c_a.Fe2N.html)
    vproto.push_back("A2B_hP12_194_cg_f;2;12;194;4;3;hP12;a,c/a,z2;C10;SiO2;beta-Tridymite");
    vproto.push_back("A4B_cI40_197_cde_c;2;20;197;4;5;cI40;a,x1,x2,x3,x4;-;Ga4Ni;Ga4Ni3");
    vproto.push_back("ABC_cP12_198_a_a_a;3;12;198;5;4;cP12;a,x1,x2,x3;F0_1;NiSSb;Ullmanite");
    vproto.push_back("A3B_cP16_198_b_a;2;16;198;4;5;cP16;a,x1,x2,y2,z2;D1;NH3;Ammonia");
    vproto.push_back("A_cP8_198_2a;1;8;198;3;3;cP8;a,x1,x2;-;N2;alpha-N2");
    vproto.push_back("AB_cP8_198_a_a;2;8;198;4;3;cP8;a,x1,x2;B21/B20/-;CO/FeSi/C1Os1;alphaCO/FeSi/C1Os1 (ICSD #168277)"); //DX20210120 - added metal-carbide
    vproto.push_back("AB_cI16_199_a_a;2;8;199;4;3;cI16;a,x1,x2;B_a;CoU;CoU");
    vproto.push_back("AB32C48_cI162_204_a_2efg_2gh;3;81;204;5;13;cI162;a,x2,x3,x4,y5,z5,y6,z6,y7,z7,x8,y8,z8;-;Mg32(Al,Zn)49;Bergman [Mg32(Al,Zn)49]");
    vproto.push_back("A3B_cI32_204_g_c;2;16;204;4;3;cI32;a,y2,z2;D0_2/-;CoAs3/O3Re1/O3Re1;Skutterudite/O3Re1 (ICSD #55462)/O3Re1 (ICSD #55469)"); //DX20210106 - added two metal-oxides
    vproto.push_back("A12B_cI26_204_g_a;2;13;204;4;3;cI26;a,y2,z2;-;Al12W;Al12W");
    vproto.push_back("A_cP8_205_c;1;8;205;3;2;cP8;a,x1;-;N2;alpha-N2");
    vproto.push_back("AB_cP16_205_c_c;2;16;205;4;3;cP16;a,x1,x2;-;CuCl;SC16");
    vproto.push_back("AB2_cP12_205_a_c;2;12;205;4;2;cP12;a,x2;C2/-;FeS2/NaO2;Pyrite/NaO2 (ICSD #87178"); //DX20200703 - added binary oxide for R. Friedrich 
    vproto.push_back("AB3C6_cI80_206_b_d_e;3;40;206;5;5;cI80;a,x2,x3,y3,z3;D5_3;(Mn,Fe)2O3;Bixbyite");
    vproto.push_back("A_cI16_206_c;1;8;206;3;2;cI16;a,x1;-;Si;BC8");
    vproto.push_back("A_cP20_213_cd;1;20;213;3;3;cP20;a,x1,y2;A13;Mn;betaMn");
    vproto.push_back("A3B4C_cP8_215_d_e_a;3;8;215;5;2;cP8;a,x3;H2_4;Cu3S4V;Sulvanite");
    vproto.push_back("AB4_cP5_215_a_e;2;5;215;4;2;cP5;a,x2;-;Fe4C;Fe4C");
    vproto.push_back("AB3C4_cP8_215_a_c_e;3;8;215;5;2;cP8;a,x3;-;Cu3AsS4;Cubic Lazarevicite");
    vproto.push_back("AB5_cF24_216_a_ce;2;6;216;4;2;cF24;a,x3;C15_b;AuBe5;AuBe5");
    vproto.push_back("ABC_cF12_216_b_c_a;3;3;216;5;1;cF12;a;C1_b;AgAsMg;Half-Heusler");
    vproto.push_back("AB_cF8_216_c_a;2;2;216;4;1;cF8;a;B3;ZnS;NaTl");
    vproto.push_back("A4B_cI10_217_c_a;2;5;217;4;2;cI10;a,x2;-;SiF4;SiF4");
    vproto.push_back("A_cI58_217_ac2g;1;29;217;3;6;cI58;a,x2,x3,z3,x4,z4;A12;Mn;alphaMn");
    vproto.push_back("A5B8_cI52_217_ce_cg;2;26;217;4;6;cI52;a,x1,x2,x3,x4,z4;-;Cu5Zn8;gammaBrass");
    vproto.push_back("A_cI16_220_c;1;8;220;3;2;cI16;a,x1;-;Li;High Pressure cI16 Li");
    vproto.push_back("A3B2_cI40_220_d_c;2;20;220;4;3;cI40;a,x1,x2;D5_c;Pu2C3;Pu2C3");
    vproto.push_back("AB_cP2_221_b_a;2;2;221;4;1;cP2;a;B2 (G0_{8});CsCl ((NH4)(NO3));CsCl and NH4NO3 I (part 3)"); //DX20210428 - added equivalent part 3 prototype (NH4NO3 I, G0_{8}, http://aflow.org/prototype-encyclopedia/AB_cP2_221_a_b.NH4.NO3.html)
    vproto.push_back("AB_cP6_221_c_d;2;6;221;4;1;cP6;a;-;NbO;NbO");
    vproto.push_back("AB3C_cP5_221_a_c_b;3;5;221;5;1;cP5;a;E2_1;CaTiO3;Cubic Perovskite");
    vproto.push_back("AB27CD3_cP32_221_a_dij_b_c;4;32;221;6;3;cP32;a,y5,y6;-;CrFe2525Ni6;Model of Austenite");
    vproto.push_back("AB3_cP4_221_a_c;2;4;221;4;1;cP4;a;L1_2;Cu3Au;Cu3Au");
    vproto.push_back("A_cP1_221_a;1;1;221;3;1;cP1;a;A_h;Po;alphaPo");
    vproto.push_back("AB11_cP36_221_c_agij;2;36;221;4;4;cP36;a,x3,y4,y5;D2_e;BaHg11;BaHg11");
    vproto.push_back("AB11CD3_cP16_221_a_dg_b_c;4;16;221;6;2;cP16;a,x5;-;CrFe11MoNi3;Model of Ferrite");
    vproto.push_back("A3B_cP4_221_d_a;2;4;221;4;1;cP4;a;D0_9;ReO3;alphaReO3");
    vproto.push_back("A6B_cP7_221_f_a;2;7;221;4;2;cP7;a,x2;D2_1;CaB6;CaB6");
    vproto.push_back("A3B_cP8_223_c_a;2;8;223;4;1;cP8;a;A15;Cr3Si;Cr3Si");
    vproto.push_back("A_cP46_223_dik;1;46;223;3;4;cP46;a,x2,y3,z3;-;Si;Si46 Clathrate");
    vproto.push_back("A2B_cP6_224_b_a;2;6;224;4;1;cP6;a;C3;Cu2O;Cuprite");
    vproto.push_back("A7B_cF32_225_bd_a;2;8;225;4;1;cF32;a;- (L1_{a});Ca7Ge (CuPt3);Ca7Ge (CuPt3)"); //DX20210428 - added equivalent part 3 prototype (CuPt3, partially occupied, http://aflow.org/prototype-encyclopedia/AB7_cF32_225_b_ad.html)
    vproto.push_back("AB3_cF16_225_a_bc;2;4;225;4;1;cF16;a;D0_3;BiF3;BiF3");
    vproto.push_back("A9B16C7_cF128_225_acd_2f_be;3;32;225;5;4;cF128;a,x5,x6,x7;-;Cr9Fe16Ni7;Model of Ferrite");
    vproto.push_back("A12B_cF52_225_i_a;2;13;225;4;2;cF52;a,y2;D2_f;UB12;UB12");
    vproto.push_back("AB2_cF12_225_a_c;2;3;225;4;1;cF12;a;C1;CaF2;Cu2Mg Cubic Laves");
    vproto.push_back("A6B23_cF116_225_e_acfh;2;29;225;4;4;cF116;a,x3,x4,y5;D8_4;Cr23C6;Cr23C6");
    vproto.push_back("AB2C_cF16_225_a_c_b;3;4;225;5;1;cF16;a;L2_1;AlCu2Mn;Heusler");
    vproto.push_back("A_cF4_225_a;1;1;225;3;1;cF4;a;A1;Cu;Face-Centered Cubic");
    vproto.push_back("AB18C8_cF108_225_a_eh_f;3;27;225;5;4;cF108;a,x2,x3,y4;-;CrFe18Ni8;Model of Austenite");
    vproto.push_back("AB_cF8_225_a_b;2;2;225;4;1;cF8;a;B1;NaCl;Rock Salt");
    vproto.push_back("A2B_cF24_227_c_a;2;6;227;4;1;cF24;a;C9;SiO2;Ideal beta-Cristobalite");
    vproto.push_back("AB2_cF96_227_e_cf;2;24;227;4;3;cF96;a,x2,x3;-;NiTi2;NiTi2");
    vproto.push_back("AB_cF16_227_a_b;2;4;227;4;1;cF16;a;B32;NaTl;NaTl");
    vproto.push_back("A_cF136_227_aeg;1;34;227;3;4;cF136;a,x2,x3,z3;-;Si;Si34 Clathrate");
    vproto.push_back("A2B_cF24_227_d_a;2;6;227;4;1;cF24;a;C15;Cu2Mg;Cu2Mg Cubic Laves");
    vproto.push_back("A_cF8_227_a;1;2;227;3;1;cF8;a;A4;C;Diamond");
    vproto.push_back("A2BC4_cF56_227_d_a_e;3;14;227;5;2;cF56;a,x3;H1_1;Al2MgO4;Spinel");
    vproto.push_back("AB2_cF48_227_c_e;2;12;227;4;2;cF48;a,x2;-;CTi2;CTi2");
    vproto.push_back("AB3C3_cF112_227_c_de_f;3;28;227;5;3;cF112;a,x3,x4;E9_3;Fe3W3C;Fe3W3C");
    vproto.push_back("A_cI2_229_a;1;1;229;3;1;cI2;a;A2;W;Body-Centered Cubic");
    vproto.push_back("A3B_cI8_229_b_a;2;4;229;4;1;cI8;a;-;H3S;High Presssure H3S");
    vproto.push_back("A4B3_cI14_229_c_b;2;7;229;4;1;cI14;a;-;Pt3O4;Pt3O4");
    vproto.push_back("A2B7_cI54_229_e_afh;2;27;229;4;4;cI54;a,x2,x3,y4;L2_2;Sb2Tl7;L22");
    vproto.push_back("AB12C3_cI32_229_a_h_b;3;16;229;5;2;cI32;a,y3;-;CrFe12Ni3;Model of Austenite");
    vproto.push_back("AB4C3_cI16_229_a_c_b;3;8;229;5;1;cI16;a;-;CrFe4Ni3;Model of Ferrite");
    vproto.push_back("A4B3_cI112_230_af_g;2;56;230;4;3;cI112;a,x2,y3;-;Ga4Ni3;Ga4Ni3");
    // -------------------------------------------------------------------------
    // Part 2
    // -------------------------------------------------------------------------
    vproto.push_back("A2B_aP6_2_aei_i;2;6;2;4;12;aP6;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3,x4,y4,z4;-;H2S;H2S");
    vproto.push_back("A8B5_mP13_6_a7b_3a2b;2;13;6;4;30;mP13;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13;-;Mo8P5;Mo8P5");
    vproto.push_back("AB_mP4_6_2b_2a;2;4;6;4;12;mP4;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;FeNi;FeNi");
    vproto.push_back("A2B_mP12_7_4a_2a;2;12;7;4;22;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;H2S;H2S IV");
    vproto.push_back("A2B_mP18_7_6a_3a;2;18;7;4;31;mP18;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;As2Ba;As2Ba");
    vproto.push_back("A3B_mP16_7_6a_2a;2;16;7;4;28;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;epsilon-WO3;epsilon-WO3");
    vproto.push_back("A9B2_mP22_7_9a_2a;2;22;7;4;37;mP22;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Rh2Ga9;Rh2Ga9");
    vproto.push_back("A5B3_mC32_9_5a_3a;2;16;9;4;28;mC32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;alpha-P3N5;alpha-P3N5");
    vproto.push_back("AB3_mC16_9_a_3a;2;8;9;4;16;mC16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;H3Cl;H3Cl");
    vproto.push_back("A2B_mP6_10_mn_bg;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;delta-PdCl2;delta-PdCl2");
    vproto.push_back("AB3_mP16_10_mn_3m3n;2;16;10;4;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;H3Cl;H3Cl");
    vproto.push_back("ABC2_mP8_10_ac_eh_mn;3;8;10;5;8;mP8;a,b/a,c/a,beta,x5,z5,x6,z6;-;AuAgTe2;Muthmannite");
    vproto.push_back("AB_mP6_10_en_am;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;LiSn;LiSn");
    vproto.push_back("A_mP8_10_2m2n;1;8;10;3;12;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;C;S-carbon");
    vproto.push_back("A7B2C2_mC22_12_aij_h_i;3;11;12;5;12;mC22;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,y5,z5;S2_{1};[Sc,Y]2Si2O7;Thortveitite");
    vproto.push_back("A_mC16_12_4i;1;8;12;3;12;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;C;M-carbon");
    vproto.push_back("A2B_mP12_13_2g_ef;2;12;13;4;12;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;H2S;H2S");
    vproto.push_back("A2B_mP6_14_e_a;2;6;14;4;7;mP6;a,b/a,c/a,beta,x2,y2,z2;-;gamma-PdCl2;gamma-PdCl2");
    vproto.push_back("A7B8_mP120_14_14e_16e;2;120;14;4;94;mP120;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30;-;alpha-C7H8;alpha-Toluene");
    vproto.push_back("AB3_mC16_15_e_cf;2;8;15;4;8;mC16;a,b/a,c/a,beta,y2,x3,y3,z3;-;H3Cl;H3Cl");
    vproto.push_back("A_mC24_15_2e2f;1;12;15;3;12;mC24;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;H;H-III");
    vproto.push_back("A2B_oP12_17_abe_e;2;12;17;4;11;oP12;a,b/a,c/a,x1,x2,x3,y3,z3,x4,y4,z4;-;alpha-Ag2Se;alpha-Naumannite");
    vproto.push_back("AB3_oP16_19_a_3a;2;16;19;4;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;H3Cl;H3Cl");
    vproto.push_back("AB2_oC6_21_a_k;2;3;21;4;4;oC6;a,b/a,c/a,z2;-/-;Ta2H/HoSb2;Ta2H/HoSb2"); //DX20210427 - added HoSb2 from part 3
    vproto.push_back("A2BC2_oF40_22_fi_ad_gh;3;10;22;5;7;oF40;a,b/a,c/a,y3,z4,z5,y6;-;CeRu2B2;CeRu2B2");
    vproto.push_back("AB_oF8_22_a_c;2;2;22;4;3;oF8;a,b/a,c/a;-;FeS;FeS");
    vproto.push_back("A3B_oI32_23_ij2k_k;2;16;23;4;14;oI32;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;H3S;H3S");
    vproto.push_back("A8B2C12D2E_oI50_23_bcfk_i_3k_j_a;5;25;23;7;18;oI50;a,b/a,c/a,x4,z5,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cu8(Fe,Zn)3Sn2S12;Stannoidite");
    vproto.push_back("ABC2_oI16_23_ab_i_k;3;8;23;5;7;oI16;a,b/a,c/a,z3,x4,y4,z4;-;NaFeS2;NaFeS2");
    vproto.push_back("ABC4_oI12_23_a_b_k;3;6;23;5;6;oI12;a,b/a,c/a,x3,y3,z3;-;BPS4;BPS4");
    vproto.push_back("AB7CD2_oI44_24_a_b3d_c_ac;4;22;24;6;17;oI44;a,b/a,c/a,x1,x2,y3,z4,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na2MgAlF7;Weberite");
    vproto.push_back("A2B_oP12_26_abc_ab;2;12;26;4;14;oP12;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5;-;H2S/beta-SeO2;H2S/beta-SeO2"); //DX20210421 - combined H2S and beta-SeO2 into one line
    vproto.push_back("A5B_oP24_26_3a3b2c_ab;2;24;26;4;25;oP24;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8,x9,y9,z9,x10,y10,z10;-;TlP5;TlP5");
    vproto.push_back("A6B4C16D_oP108_27_abcd4e_4e_16e_e;4;108;27;6;82;oP108;a,b/a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29;-;Ca4Al6O16S;Ca4Al6O16S");
    vproto.push_back("A2B_oP12_29_2a_a;2;12;29;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-/-;ZrO2/O2Zr1;ZrO2/O2Zr1 (ICSD #67004)"); //DX20210106 - added metal-oxide prototype
    vproto.push_back("AB2_oP12_29_a_2a;2;12;29;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;FeS2;Pyrite");
    vproto.push_back("ABC_oP12_29_a_a_a;3;12;29;5;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;CoAsS;Cobaltite");
    vproto.push_back("A5B3C15_oP46_30_a2c_bc_a7c;3;46;30;5;36;oP46;a,b/a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;Bi5Nb3O15;Bi5Nb3O15");
    vproto.push_back("ABC3_oP20_30_2a_c_3c;3;20;30;5;17;oP20;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CuBrSe3;CuBrSe3");
    vproto.push_back("A13B2C2_oP34_32_a6c_c_c;3;34;32;5;28;oP34;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Re2O5[SO4]2;Re2O52");
    vproto.push_back("A2B3_oP40_33_4a_6a;2;40;33;4;33;oP40;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;kappa-Al2O3;kappa-alumina");
    vproto.push_back("A2B8C_oP22_34_c_4c_a;3;22;34;5;19;oP22;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;TiAl2Br8;TiAl2Br8");
    vproto.push_back("AB2_oP6_34_a_c;2;6;34;4;7;oP6;a,b/a,c/a,z1,x2,y2,z2;-;FeSb2;FeSb2");
    vproto.push_back("AB8C2_oC22_35_a_ab3e_e;3;11;35;5;14;oC22;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6,y7,z7;-/-;V2MoO8/Mo1O8V2;V2MoO8/Mo1O8V2 (ICSD #25378)");
    vproto.push_back("AB_oC8_36_a_a;2;4;36;4;7;oC8;a,b/a,c/a,y1,z1,y2,z2;-;HCl;HCl");
    vproto.push_back("A2B5C2_oC36_37_d_c2d_d;3;18;37;5;16;oC36;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Li2Si2O5;Li2Si2O5");
    vproto.push_back("A2B3_oC40_39_2d_2c2d;2;20;39;4;19;oC40;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Ta3S2;Ta3S2");
    vproto.push_back("A9BC_oC44_39_3c3d_a_c;3;22;39;5;21;oC44;a,b/a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;VPCl9;VPCl9");
    vproto.push_back("AB2C_oC16_40_a_2b_b;3;8;40;5;10;oC16;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-;K2CdPb;K2CdPb");
    vproto.push_back("AB3_oC16_40_b_3b;2;8;40;4;11;oC16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-;CeTe3;CeTe3");
    vproto.push_back("A10B3_oF52_42_2abce_ab;2;13;42;4;13;oF52;a,b/a,c/a,z1,z2,z3,z4,z5,y6,z6,x7,y7,z7;-;W3O10;W3O10");
    vproto.push_back("AB_oF8_42_a_a;2;2;42;4;5;oF8;a,b/a,c/a,z1,z2;-;BN;BN");
    vproto.push_back("A2BC2_oI20_45_c_b_c;3;10;45;5;10;oI20;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;-;MnGa2Sb2;MnGa2Sb2");
    vproto.push_back("ABC_oI36_46_ac_bc_3b;3;18;46;5;18;oI36;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4,y5,z5,x6,y6,z6,x7,y7,z7;-;TiFeSi;TiFeSi");
    vproto.push_back("A2B8CD_oP24_48_k_2m_d_b;4;24;48;6;10;oP24;a,b/a,c/a,z3,x4,y4,z4,x5,y5,z5;-;alpha-RbPr[MoO4]2;alpha-RbPr2");
    vproto.push_back("A5B2_oP14_49_dehq_ab;2;14;49;4;5;oP14;a,b/a,c/a,x6,y6;-;beta-Ta2O5;beta-Ta2O5");
    vproto.push_back("AB2C8D_oP24_49_g_q_2qr_e;4;24;49;6;12;oP24;a,b/a,c/a,x3,y3,x4,y4,x5,y5,x6,y6,z6;-;CsPr[MoO4]2;CsPr2");
    vproto.push_back("A2BC4_oP28_50_ij_ac_ijm;3;28;50;5;10;oP28;a,b/a,c/a,y3,y4,y5,y6,x7,y7,z7;-;La2NiO4;La2NiO4");
    vproto.push_back("A3BC2_oP48_50_3m_m_2m;3;48;50;5;21;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-Tl2TeO3;alpha-Tl2TeO3");
    vproto.push_back("A2B_oP24_52_2e_cd;2;24;52;4;11;oP24;a,b/a,c/a,z1,x2,x3,y3,z3,x4,y4,z4;-;GaCl2;GaCl2");
    vproto.push_back("A3B2_oP20_52_de_cd;2;20;52;4;9;oP20;a,b/a,c/a,z1,x2,x3,x4,y4,z4;-;Sr2Bi3;Sr2Bi3");
    vproto.push_back("ABC2_oP16_53_h_e_gh;3;16;53;5;9;oP16;a,b/a,c/a,x1,y2,y3,z3,y4,z4;-;TaNiTe2;TaNiTe2");
    vproto.push_back("ABC3_oP20_53_e_g_hi;3;20;53;5;10;oP20;a,b/a,c/a,x1,y2,y3,z3,x4,y4,z4;-;CuBrSe3;CuBrSe3");
    vproto.push_back("ABC3_oP20_54_e_d_cf;3;20;54;5;9;oP20;a,b/a,c/a,y1,z2,z3,x4,y4,z4;-;BiGaO3;BiGaO3");
    vproto.push_back("A2B_oP24_55_2g2h_gh;2;24;55;4;15;oP24;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6;-;GeAs2;GeAs2");
    vproto.push_back("A3B5_oP16_55_ch_agh;2;16;55;4;9;oP16;a,b/a,c/a,x3,y3,x4,y4,x5,y5;-;Rh5Ge3;Rh5Ge3");
    vproto.push_back("A_oP16_55_2g2h;1;16;55;3;11;oP16;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4;-;C;R-carbon");
    vproto.push_back("A2B_oP6_58_g_a;2;6;58;4;5;oP6;a,b/a,c/a,x2,y2;C50/-/-/-/-;alpha-PdCl2/N2Os1/N2Pt1/Ir1N2/Au1N2;alpha-PdCl2/N2Os1 (ICSD #157283)/N2Pt1 (ICSD #166463)/Ir1N2 (ICSD #160620)/Au1N2 (ICSD #166465)"); //DX20201019 - added four metal nitrides
    vproto.push_back("ABC_oP6_59_a_b_a;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;E0_{5};FeOCl;FeOCl"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A2B3_oP20_60_d_cd;2;20;60;4;10;oP20;a,b/a,c/a,y1,x2,y2,z2,x3,y3,z3;-;Rh2S3;Rh2S3");
    vproto.push_back("A3B_oP32_60_3d_d;2;32;60;4;15;oP32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;WO3;WO3");
    vproto.push_back("A7B8_oP120_60_7d_8d;2;120;60;4;48;oP120;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;beta-C7H8;beta-Toluene");
    vproto.push_back("AB_oP48_61_3c_3c;2;48;61;4;21;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Benzene;Benzene");
    vproto.push_back("A2B3_oP20_62_2c_3c;2;20;62;4;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;D5_{10};Cr3C2;Tongbaite");
    vproto.push_back("A2B4C_oP28_62_ac_2cd_c;3;28;62;5;14;oP28;a,b/a,c/a,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;S1_{2}/-;Mg2SiO4/Ag2CrO4;Forsterite/Ag2CrO4 (ICSD #16298)"); //DX20200703 - added R. Friedrich ternary oxide
    //DX 20180619 - added info to label in part 1 vproto.push_back("A2B_oP12_62_2c_c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;C29;SrH2;SrH2");
    vproto.push_back("A3B_oP16_62_cd_c;2;16;62;4;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;D0_{20}/-;epsilon-NiAl3/Fe3N1;epsilon-NiAl3/Fe3N1 (ICSD #260758)"); //DX20201103 - added metal-nitride
    vproto.push_back("AB2C3_oP24_62_c_d_cd;3;24;62;5;13;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3,x4,y4,z4;E9_{e};CuFe2S3;Cubanite");
    vproto.push_back("AB3_oP16_62_c_3c;2;16;62;4;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;D0_{8}/-;MoO3/Mo1O3;Molybdite/Mo1O3 (ICSD #644065"); //DX20210106 - added metal-nitride
    vproto.push_back("AB4C_oP24_62_c_2cd_c;3;24;62;5;14;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5;H0_{2};BaSO4;Barite");
    //DX 20180619 - added info to label in part 1 vproto.push_back("AB_oP8_62_c_c;2;8;62;4;7;oP8;a,b/a,c/a,x1,z1,x2,z2;B14;FeAs;Westerveldite");
    vproto.push_back("A2BC3_oC24_63_e_c_cg;3;12;63;5;8;oC24;a,b/a,c/a,y1,y2,x3,x4,y4;-;KFe2S3;Rasvumite");
    vproto.push_back("A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h;3;130;63;5;59;oC260;a,b/a,c/a,y1,y2,y3,x4,y5,z5,y6,z6,y7,z7,y8,z8,y9,z9,y10,z10,y11,z11,y12,z12,y13,z13,y14,z14,y15,z15,y16,z16,x17,y17,x18,y18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26;-;La43Ni17Mg5;La43Ni17Mg5");
    vproto.push_back("A6B_oC28_63_efg_c;2;14;63;4;9;oC28;a,b/a,c/a,y1,x2,y3,z3,x4,y4;D2_{h};MnAl6;MnAl6");
    vproto.push_back("AB3C_oC20_63_a_cf_c;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;MgSiO3 (AsCV3);Post-perovskite and V3AsC (part 3)"); //DX20210428 - added equivalent part 3 prototype (V3AsC, http://aflow.org/prototype-encyclopedia/ABC3_oC20_63_c_b_cf.html, part 3)
    vproto.push_back("AB4C_oC24_63_a_fg_c;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,z3,x4,y4;-/-/-;MgO4S/In1O4V1/Cr1O4V1;MgSO4/In1O4V1 (ICSD #10431)/Cr1O4V1 (ICSD #27508)"); //DX20210106 - added two metal-oxide
    vproto.push_back("AB4C_oC24_63_c_fg_c;3;12;63;5;9;oC24;a,b/a,c/a,y1,y2,y3,z3,x4,y4;H0_{1};CaSO4;Anhydrite");
    vproto.push_back("A2B_oC24_64_2f_f;2;12;64;4;9;oC24;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-;H2S;H2S");
    vproto.push_back("A2B4C_oC28_66_l_kl_a;3;14;66;5;8;oC28;a,b/a,c/a,z2,x3,y3,x4,y4;-;SrAl2Se4;SrAl2Se4");
    vproto.push_back("A3B_oC64_66_gi2lm_2l;2;32;66;4;16;oC64;a,b/a,c/a,x1,z2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;-;H3S;H3S");
    vproto.push_back("A3B_oC64_66_kl2m_bdl;2;32;66;4;14;oC64;a,b/a,c/a,z3,x4,y4,x5,y5,x6,y6,z6,x7,y7,z7;-;beta-ThI3;beta-ThI3");
    vproto.push_back("A2BC_oC16_67_ag_b_g;3;8;67;5;5;oC16;a,b/a,c/a,z3,z4;-;Al2CuIr;Al2CuIr");
    vproto.push_back("ABC2_oC16_67_b_g_ag;3;8;67;5;5;oC16;a,b/a,c/a,z3,z4;-;HoCuP2;HoCuP2");
    vproto.push_back("AB_oC8_67_a_g;2;4;67;4;4;oC8;a,b/a,c/a,z2;-;alpha-FeSe/alpha-PbO;alpha-FeSe/alpha-PbO"); //DX20210421 - combined alpha-FeSe/PbO into one line
    vproto.push_back("AB4_oC20_68_a_i;2;10;68;4;6;oC20;a,b/a,c/a,x2,y2,z2;-;PdSn4;PdSn4");
    vproto.push_back("AB2_oF48_70_f_fg;2;12;70;4;6;oF48;a,b/a,c/a,y1,y2,z3;D1_{f} (C_{b});Mn2B (Mg2Cu);Mn2B (Mg2Cu)"); //DX20210428 - added similar part 3 prototype with different Struk. designation (Mg2Cu, C_{b}, http://aflow.org/prototype-encyclopedia/AB2_oF48_70_g_fg.html)
    vproto.push_back("A4B3_oI14_71_gh_cg;2;7;71;4;6;oI14;a,b/a,c/a,y2,y3,y4;D7_{b};Ta3B4;Ta3B4");
    vproto.push_back("ABC_oI12_71_h_j_g;3;6;71;5;6;oI12;a,b/a,c/a,y1,y2,z3;-;NbPS;NbPS");
    vproto.push_back("ABCD3_oI48_73_d_e_e_ef;4;24;73;6;10;oI48;a,b/a,c/a,y1,z2,z3,z4,x5,y5,z5;-;KAg[CO3];KAg");
    vproto.push_back("A2B_oI12_74_h_e;2;6;74;4;6;oI12;a,b/a,c/a,z1,y2,z2;-;KHg2 (CeCu2);KHg2 (CeCu2)"); //DX20210428 - added equivalent part 3 prototype (CeCu2, http://aflow.org/prototype-encyclopedia/AB2_oI12_74_e_h.html)
    vproto.push_back("A4B_oI20_74_beh_e;2;10;74;4;7;oI20;a,b/a,c/a,z2,z3,y4,z4;D1_{b};Al4U;Al4U");
    vproto.push_back("AB2C12D4_tP76_75_2a2b_2d_12d_4d;4;76;75;6;60;tP76;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;BaCr2Ru4O12;BaCr2Ru4O12");
    vproto.push_back("A2BC_tP16_76_2a_a_a;3;16;76;5;14;tP16;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-;LaRhC2/C2Ce1Rh1;LaRhC2/C2Ce1Rh1 (ICSD #67327)"); //DX20210121 - added metal-carbide
    vproto.push_back("A3B7_tP40_76_3a_7a;2;40;76;4;32;tP40;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cs3P7;Cs3P7");
    vproto.push_back("A2B6CD7_tP64_77_2d_6d_d_ab6d;4;64;77;6;49;tP64;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;-;MgB2O(OH)6;Pinnoite");
    vproto.push_back("A2B_tP48_77_8d_4d;2;48;77;4;38;tP48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;H2S;H2S III");
    vproto.push_back("A2B7C2_tP88_78_4a_14a_4a;3;88;78;5;68;tP88;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;Sr2As2O7;Sr2As2O7");
    vproto.push_back("A2BC2_tI20_79_c_2a_c;3;10;79;5;10;tI20;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4;-;TlZn2Sb2;TlZn2Sb2");
    vproto.push_back("AB2_tI48_80_2b_4b;2;24;80;4;20;tI48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;beta-NbO2;beta-NbO2");
    vproto.push_back("AB2_tP12_81_adg_2h;2;12;81;4;9;tP12;a,c/a,z3,x4,y4,z4,x5,y5,z5;-;GeSe2;GeSe2");
    vproto.push_back("A3B_tI32_82_3g_g;2;16;82;4;14;tI32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;D0_{e};Ni3P;Ni3P");
    vproto.push_back("A3B2_tP10_83_adk_j;2;10;83;4;6;tP10;a,c/a,x3,y3,x4,y4;-;Ti2Ge3;Ti2Ge3");
    vproto.push_back("A2B_tP30_85_ab2g_cg;2;30;85;4;12;tP30;a,c/a,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;SrBr2;SrBr2");
    vproto.push_back("AB3_tP32_86_g_3g;2;32;86;4;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Ti3P;Ti3P");
    vproto.push_back("A4B_tI20_88_f_a;2;10;88;4;5;tI20;a,c/a,x2,y2,z2;-;ThCl4;ThCl4");
    vproto.push_back("AB2_tI96_88_2f_4f;2;48;88;4;20;tI96;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-NbO2;alpha-NbO2");
    vproto.push_back("A17BC4D_tP184_89_17p_p_4p_io;4;184;89;6;70;tP184;a,c/a,z1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24;-;C17FeO4Pt;C17FeO4Pt");
    vproto.push_back("A4B2C13D_tP40_90_g_d_cef2g_c;4;40;90;6;16;tP40;a,c/a,z1,z2,z3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na4Ti2Si8O22[H2O]4;Na4Ti2Si8O224");
    vproto.push_back("AB4C17D4E_tP54_90_a_g_c4g_g_c;5;54;90;7;22;tP54;a,c/a,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;BaCu4[VO][PO4]4;BaCu44");
    vproto.push_back("ABC_tP24_91_d_d_d;3;24;91;5;11;tP24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ThBC;ThBC");
    vproto.push_back("AB32CD4E8_tP184_93_i_16p_af_2p_4p;5;184;93;7;69;tP184;a,c/a,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25;-;AsPh4CeS8P4Me8;AsPh4CeS8P4Me8");
    vproto.push_back("A14B3C5_tP44_94_c3g_ad_bg;3;44;94;5;16;tP44;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na5Fe3F14;Na5Fe3F14");
    vproto.push_back("A6B2C_tP18_94_eg_c_a;3;18;94;5;7;tP18;a,c/a,z2,x3,x4,y4,z4;-;Li2MoF6;Li2MoF6");
    vproto.push_back("ABC_tP24_95_d_d_d;3;24;95;5;11;tP24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;ThBC;ThBC");
    vproto.push_back("A2B8CD_tI24_97_d_k_a_b;4;12;97;6;5;tI24;a,c/a,x4,y4,z4;-;NaGdCu2F8;NaGdCu2F8");
    vproto.push_back("AB8C2_tI44_97_e_2k_cd;3;22;97;5;9;tI44;a,c/a,z3,x4,y4,z4,x5,y5,z5;-;Ta2Se8I;Ta2Se8I");
    vproto.push_back("A2B_tI12_98_f_a;2;6;98;4;3;tI12;a,c/a,x2;-;CdAs2;CdAs2");
    vproto.push_back("A2B8C2D_tP26_100_c_abcd_c_a;4;26;100;6;14;tP26;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;Ba2TiSi2O8;Fresnoite");
    vproto.push_back("A3B11C6_tP40_100_ac_bc2d_cd;3;40;100;5;19;tP40;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Ce3Si6N11;Ce3Si6N11");
    vproto.push_back("A7B7C2_tP32_101_bde_ade_d;3;32;101;5;16;tP32;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;gamma-MgNiSn;gamma-MgNiSn");
    vproto.push_back("A2B3_tP20_102_2c_b2c;2;20;102;4;11;tP20;a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Gd3Al2;Gd3Al2");
    vproto.push_back("AB4_tP10_103_a_d;2;10;103;4;6;tP10;a,c/a,z1,x2,y2,z2;-;NbTe4;NbTe4");
    vproto.push_back("A5B5C4_tP28_104_ac_ac_c;3;28;104;5;13;tP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Ba5In4Bi5;Ba5In4Bi5");
    vproto.push_back("AB6C4_tP22_104_a_2ac_c;3;22;104;5;11;tP22;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5;-;Tl4HgI6;Tl4HgI6");
    vproto.push_back("A2BC2_tP20_105_f_ac_2e;3;20;105;5;11;tP20;a,c/a,z1,z2,x3,z3,x4,z4,x5,y5,z5;-;BaGe2As2;BaGe2As2");
    vproto.push_back("A3BC3D_tP64_106_3c_c_3c_c;4;64;106;6;26;tP64;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;NaZn[OH]3;NaZn3");
    vproto.push_back("A5B7_tI24_107_ac_abd;2;12;107;4;9;tI24;a,c/a,z1,z2,z3,x4,z4,x5,z5;-;Co5Ge7;Co5Ge7");
    vproto.push_back("AB_tI4_107_a_a;2;2;107;4;4;tI4;a,c/a,z1,z2;-;GeP;GeP");
    vproto.push_back("A3B5_tI32_108_ac_a2c;2;16;108;4;10;tI32;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Sr5Si3;Sr5Si3");
    vproto.push_back("ABC_tI12_109_a_a_a;3;6;109;5;5;tI12;a,c/a,z1,z2,z3;-;LaPtSi;LaPtSi");
    vproto.push_back("AB_tI8_109_a_a;2;4;109;4;4;tI8;a,c/a,z1,z2;-;NbAs;NbAs");
    vproto.push_back("A2BC8_tI176_110_2b_b_8b;3;88;110;5;35;tI176;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Be[BH4]2;Be2");
    vproto.push_back("A2B_tP12_111_2n_adf;2;12;111;4;6;tP12;a,c/a,x4,z4,x5,z5;-;MnF2;MnF2");
    vproto.push_back("AB_tP8_111_n_n;2;8;111;4;6;tP8;a,c/a,x1,z1,x2,z2;-;VN;VN"); //DX 20180925 - prototype name should be VN not NV
    vproto.push_back("AB4C_tP12_112_b_n_e;3;12;112;5;5;tP12;a,c/a,x3,y3,z3;-;alpha-CuAlCl4;alpha-CuAlCl4");
    vproto.push_back("A2BC7D2_tP24_113_e_a_cef_e;4;24;113;6;12;tP24;a,c/a,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;S5_{3};Ca2MgSi2O7;Akermanite");
    vproto.push_back("A3B_tP32_114_3e_e;2;32;114;4;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;SeO3;SeO3");
    vproto.push_back("A4B_tP10_114_e_a;2;10;114;4;5;tP10;a,c/a,x2,y2,z2;-;Pd4Se;Pd4Se");
    vproto.push_back("A2B3_tP5_115_g_ag;2;5;115;4;4;tP5;a,c/a,z2,z3;-;Rh3P2;Rh3P2");
    vproto.push_back("AB2_tP12_115_j_egi;2;12;115;4;7;tP12;a,c/a,z1,z2,x3,x4,z4;-;HgI2;HgI2");
    vproto.push_back("A2B3_tP20_116_bci_fj;2;20;116;4;7;tP20;a,c/a,x3,z4,x5,y5,z5;-;Ru2Sn3;Ru2Sn3");
    vproto.push_back("A2B3_tP20_117_i_adgh;2;20;117;4;7;tP20;a,c/a,x3,x4,x5,y5,z5;D5_{12};beta-Bi2O3;beta-Bi2O3"); //DX 20180925 - added Strukturbericht designation
    vproto.push_back("A3B_tP16_118_ei_f;2;16;118;4;7;tP16;a,c/a,z1,x2,x3,y3,z3;-;RuIn3;RuIn3");
    vproto.push_back("A5B3_tP32_118_g2i_aceh;2;32;118;4;11;tP32;a,c/a,z3,x4,z5,x6,y6,z6,x7,y7,z7;-;Ir3Ga5;Ir3Ga5");
    vproto.push_back("A3B_tI24_119_b2i_af;2;12;119;4;7;tI24;a,c/a,z3,x4,z4,x5,z5;-;RbGa3;RbGa3");
    vproto.push_back("AB_tI4_119_c_a;2;2;119;4;2;tI4;a,c/a;-;GaSb;GaSb");
    vproto.push_back("A4BC2_tI28_120_i_d_e;3;14;120;5;6;tI28;a,c/a,x2,x3,y3,z3;-;KAu4Sn2;KAu4Sn2");
    vproto.push_back("A4BC4D_tP10_123_gh_a_i_d;4;10;123;6;5;tP10;a,c/a,z3,z4,z5;-;CaRbFe4As4;CaRbFe4As4"); //DX 20180925 - prototype name should be Ca not Cs
    vproto.push_back("AB4C_tP12_124_a_m_c;3;12;124;5;4;tP12;a,c/a,x3,y3;-;Nb4CoSi;Nb4CoSi");
    vproto.push_back("AB4_tP10_124_a_m;2;10;124;4;4;tP10;a,c/a,x2,y2;-;NbTe4;NbTe4");
    vproto.push_back("A4B_tP10_125_m_a;2;10;125;4;4;tP10;a,c/a,x2,z2;-;PtPb4;PtPb4");
    vproto.push_back("ABC4_tP12_125_a_b_m;3;12;125;5;4;tP12;a,c/a,x3,z3;-;KCeSe4;KCeSe4");
    vproto.push_back("A2BC4_tP28_126_cd_e_k;3;28;126;5;6;tP28;a,c/a,z3,x4,y4,z4;-;BiAl2S4;BiAl2S4");
    vproto.push_back("A4B_tP20_127_ehj_g;2;20;127;4;7;tP20;a,c/a,z1,x2,x3,x4,y4;D1_{e};ThB4;ThB4");
    vproto.push_back("A6B2C_tP18_128_eh_d_a;3;18;128;5;5;tP18;a,c/a,z3,x4,y4;-;K2SnCl6;K2SnCl6");
    vproto.push_back("A7B2C_tP40_128_egi_h_e;3;40;128;5;10;tP40;a,c/a,z1,z2,x3,x4,y4,x5,y5,z5;E9_{a};FeCu2Al7;FeCu2Al7");
    vproto.push_back("A2BC4_tP28_130_f_c_g;3;28;130;5;7;tP28;a,c/a,z1,x2,x3,y3,z3;-;CuBi2O4;CuBi2O4");
    vproto.push_back("A5B3_tP32_130_cg_cf;2;32;130;4;8;tP32;a,c/a,z1,z2,x3,x4,y4,z4;-;Ba5Si3;Ba5Si3");
    vproto.push_back("A2B2C4D_tP18_132_e_i_o_d;4;18;132;6;5;tP18;a,c/a,x3,x4,z4;-;Rb2TiCu2Se4;Rb2TiCu2S4");
    vproto.push_back("AB6C_tP16_132_d_io_a;3;16;132;5;5;tP16;a,c/a,x3,x4,z4;-;AgUF6;AgUF6");
    vproto.push_back("AB3_tP32_133_h_i2j;2;32;133;4;6;tP32;a,c/a,x1,x2,x3,x4;-;beta-V3S;beta-V3S");
    vproto.push_back("A2B_tP24_135_gh_h;2;24;135;4;7;tP24;a,c/a,x1,x2,y2,x3,y3;C47;SeO2;Downeyite");
    vproto.push_back("A4B2C_tP28_135_gh_h_d;3;28;135;5;7;tP28;a,c/a,x2,x3,y3,x4,y4;-;ZnSb2O4;ZnSb2O4");
    vproto.push_back("A2B3_tP40_137_cdf_3g;2;40;137;4;11;tP40;a,c/a,z1,z2,x3,y4,z4,y5,z5,y6,z6;D5_{9};Zn3P2;Zn3P2");
    vproto.push_back("A2B_tP6_137_d_a;2;6;137;4;3;tP6;a,c/a,z2;-/-;ZrO2/O2Zr1;ZrO2/O2Zr1 (ICSD #93028)"); //DX20210106 - added metal-oxide prototype
    vproto.push_back("A4BC4_tP18_137_g_b_g;3;18;137;5;6;tP18;a,c/a,y2,z2,y3,z3;-;CeCo4B4;CeCo4B4");
    vproto.push_back("AB2_tP6_137_a_d;2;6;137;4;3;tP6;a,c/a,z2;C13;HgI2;HgI2");
    vproto.push_back("A_tP12_138_bi;1;12;138;3;4;tP12;a,c/a,x2,z2;-;C;C");
    vproto.push_back("AB_tI8_139_e_e;2;4;139;4;4;tI8;a,c/a,z1,z2;D3_{1}/-;Hg2Cl2/-;Calomel/-");
    vproto.push_back("A3B5_tI32_140_ah_bk;2;16;140;4;5;tI32;a,c/a,x3,x4,y4;D8_{m};W5Si3;W5Si3");
    vproto.push_back("A3B5_tI32_140_ah_cl;2;16;140;4;5;tI32;a,c/a,x3,x4,z4;D8_{l};Cr5B3;Cr5B3");
    //DX 20180619 - added info to label in part 1 vproto.push_back("A2B_tI12_141_e_a;2;6;141;4;3;tI12;a,c/a,z2;C_{c};alpha-ThSi2;alpha-ThSi2");
    vproto.push_back("A_tI16_142_f;1;8;142;3;3;tI16;a,c/a,x1;-;S;S-III");
    vproto.push_back("A4B14C3_hP21_143_bd_ac4d_d;3;21;143;5;23;hP21;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Ta3Al4O13[OH];Simpsonite");
    vproto.push_back("A4B6C_hP11_143_bd_2d_a;3;11;143;5;13;hP11;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;ScRh6P4;ScRh6P4");
    vproto.push_back("AB2_hP12_143_cd_ab2d;2;12;143;4;14;hP12;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;MoS2;MoS2");
    vproto.push_back("A4B_hP15_144_4a_a;2;15;144;4;17;hP15;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;IrGe4;IrGe4");
    vproto.push_back("AB_hP6_144_a_a;2;6;144;4;8;hP6;a,c/a,x1,y1,z1,x2,y2,z2;-;ZnTe;ZnTe"); //DX 20180925 - prototype name should be ZnTe not TeZn
    vproto.push_back("A2B3C3DE7_hP48_145_2a_3a_3a_a_7a;5;48;145;7;50;hP48;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;NaCa3[CO3]2F3[H2O];Sheldrickite");
    vproto.push_back("A3BC_hR5_146_b_a_a;3;5;146;5;7;hR5;a,c/a,x1,x2,x3,y3,z3;-;gamma-Ag3SI;gamma-Ag3SI");
    vproto.push_back("ABC3_hR10_146_2a_2a_2b;3;10;146;5;12;hR10;a,c/a,x1,x2,x3,x4,x5,y5,z5,x6,y6,z6;-;FePSe3;FePSe3");
    vproto.push_back("A2B4C_hR42_148_2f_4f_f;3;42;148;5;23;hR42;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;S1_{3};Be2SiO4;Phenakite");
    vproto.push_back("A2B_hR18_148_2f_f;2;18;148;4;11;hR18;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;beta-PdCl2;beta-PdCl2");
    vproto.push_back("AB3_hP24_149_acgi_3l;2;24;149;4;13;hP24;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Ti3O;Ti3O");
    vproto.push_back("A3B_hP24_153_3c_2b;2;24;153;4;13;hP24;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;CrCl3;CrCl3");
    vproto.push_back("A_hP9_154_bc;1;9;154;3;6;hP9;a,c/a,x1,x2,y2,z2;-;S;S-II");
    vproto.push_back("AB2_hP9_156_b2c_3a2bc;2;9;156;4;11;hP9;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;-;CdI2;CdI2");
    vproto.push_back("AB_hP12_156_2ab3c_2ab3c;2;12;156;4;14;hP12;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12;-;CuI;CuI");
    vproto.push_back("AB_hP4_156_ab_ab;2;4;156;4;6;hP4;a,c/a,z1,z2,z3,z4;-;beta-CuI;beta-CuI");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB_hP4_156_ac_ac;2;4;156;4;6;hP4;a,c/a,z1,z2,z3,z4;-;beta-CuI;beta-CuI");
    vproto.push_back("A5B6C2_hP13_157_2ac_2c_b;3;13;157;5;11;hP13;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,z6;-;Ag5Pb2O6;Ag5Pb2O6");
    vproto.push_back("A3B_hP8_158_d_a;2;8;158;4;6;hP8;a,c/a,z1,x2,y2,z2;-;beta-RuCl3;beta-RuCl3");
    vproto.push_back("A2B3_hP20_159_bc_2c;2;20;159;4;12;hP20;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Bi2O3;Bi2O3");
    vproto.push_back("A4B3_hP28_159_ab2c_2c;2;28;159;4;16;hP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;alpha-Si3N4;Nierite");
    vproto.push_back("AB4C7D_hP26_159_b_ac_a2c_b;4;26;159;6;15;hP26;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;YbBaCo4O7;YbBaCo4O7");
    vproto.push_back("A3B_hR4_160_b_a;2;4;160;4;5;hR4;a,c/a,x1,x2,z2;-;H3S;H3S");
    vproto.push_back("A8B5_hR26_160_a3bc_a3b;2;26;160;4;19;hR26;a,c/a,x1,x2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9;D8_{10};Cr5Al8;Cr5Al8"); //DX 20180925 - prototype name should be Cr5Al8 not Al8Cr5
    vproto.push_back("ABC_hR3_160_a_a_a;3;3;160;5;5;hR3;a,c/a,x1,x2,x3;F0_{2}/-;COS/Ag1C1N1;Carbonyl Sulphide/Ag1C1N1 (ICSD #85783)"); //DX20210221 - added carbo-nitride
    vproto.push_back("AB_hR10_160_5a_5a;2;10;160;4;12;hR10;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;B7;SiC;Moissanite-15R");
    vproto.push_back("A2B3_hP5_164_d_ad;2;5;164;4;4;hP5;a,c/a,z2,z3;D5_{2};La2O3;La2O3");
    vproto.push_back("AB2_hP9_164_bd_c2d;2;9;164;4;6;hP9;a,c/a,z2,z3,z4,z5;-;deltaH^II-NW2;deltaH^II-NW2");
    vproto.push_back("ABC2_hP4_164_a_b_d;3;4;164;5;3;hP4;a,c/a,z3;-;CuNiSb2;CuNiSb2");
    vproto.push_back("A3B_hP24_165_bdg_f;2;24;165;4;7;hP24;a,c/a,z2,x3,x4,y4,z4;D0_{21};Cu3P;Cu3P");
    vproto.push_back("A4B3_hR7_166_2c_ac;2;7;166;4;5;hR7;a,c/a,x2,x3,x4;D7_{1};Al4C3;Al4C3");
    vproto.push_back("ABC_hR6_166_c_c_c;3;6;166;5;5;hR6;a,c/a,x1,x2,x3;-;SmSI;SmSI");
    vproto.push_back("AB3C_hR10_167_b_e_a;3;10;167;5;3;hR10;a,c/a,x3;-;PrNiO3;PrNiO3");
    vproto.push_back("ABC2_hR24_167_e_e_2e;3;24;167;5;6;hR24;a,c/a,x1,x2,x3,x4;F5_{13};KBO2;KBO2");
    vproto.push_back("A2B13C4_hP57_168_d_c6d_2d;3;57;168;5;30;hP57;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;K2Ta4O9F4;K2Ta4O9F4");
    vproto.push_back("AB4C_hP72_168_2d_8d_2d;3;72;168;5;38;hP72;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Al[PO4];Al");
    vproto.push_back("A2B3_hP30_169_2a_3a;2;30;169;4;17;hP30;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;alpha-Al2S3;alpha-Al2S3");
    vproto.push_back("A2B3_hP30_170_2a_3a;2;30;170;4;17;hP30;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Al2S3;Al2S3");
    vproto.push_back("A10B2C_hP39_171_5c_c_a;3;39;171;5;21;hP39;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Sr[S2O6][H2O]4;Sr4");
    vproto.push_back("A10B2C_hP39_172_5c_c_a;3;39;172;5;21;hP39;a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Sr[S2O6][H2O]4;Sr4");
    vproto.push_back("A3B_hP8_173_c_b;2;8;173;4;6;hP8;a,c/a,z1,x2,y2,z2;-;PI3;PI3");
    vproto.push_back("A4B3_hP14_173_bc_c;2;14;173;4;9;hP14;a,c/a,z1,x2,y2,z2,x3,y3,z3;-;beta-Si3N4;beta-Si3N4");
    vproto.push_back("A12B7C2_hP21_174_2j2k_ajk_cf;3;21;174;5;14;hP21;a,c/a,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9;-;Fe12Zr2P7;Fe12Zr2P7");
    vproto.push_back("ABC_hP12_174_cj_fk_aj;3;12;174;5;8;hP12;a,c/a,x4,y4,x5,y5,x6,y6;-;GdSI;GdSI");
    vproto.push_back("A8B7C6_hP21_175_ck_aj_k;3;21;175;5;8;hP21;a,c/a,x3,y3,x4,y4,x5,y5;-;Nb7Ru6B8;Nb7Ru6B8");
    vproto.push_back("ABC_hP36_175_jk_jk_jk;3;36;175;5;14;hP36;a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6;-;Mg[NH];Mg");
    vproto.push_back("A3B2_hP10_176_h_bc;2;10;176;4;4;hP10;a,c/a,x3,y3;-;Er3Ru2;Er3Ru2");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B2_hP10_176_h_bd;2;10;176;4;4;hP10;a,c/a,x3,y3;-;Er3Ru2;Er3Ru2");
    vproto.push_back("A3B3C_hP14_176_h_h_c;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Fe3Te3Tl;Fe3Te3Tl");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B3C_hP14_176_h_h_d;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Fe3Te3Tl;Fe3Te3Tl");
    vproto.push_back("A3B_hP8_176_h_c;2;8;176;4;4;hP8;a,c/a,x2,y2;-;UCl3;UCl3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3B_hP8_176_h_d;2;8;176;4;4;hP8;a,c/a,x2,y2;-;UCl3;UCl3");
    vproto.push_back("A2B_hP36_177_j2lm_n;2;36;177;4;9;hP36;a,c/a,x1,x2,x3,x4,x5,y5,z5;-;SiO2;SiO2");
    vproto.push_back("AB3_hP24_178_b_ac;2;24;178;4;7;hP24;a,c/a,x1,x2,x3,y3,z3;-;AuF3;AuF3");
    vproto.push_back("A_hP6_178_a;1;6;178;3;3;hP6;a,c/a,x1;-;Sc;Sc-V");
    vproto.push_back("AB3_hP24_179_b_ac;2;24;179;4;7;hP24;a,c/a,x1,x2,x3,y3,z3;-;AuF3;AuF3");
    vproto.push_back("A2B_hP9_181_j_c;2;9;181;4;3;hP9;a,c/a,x2;-;beta-SiO2;beta-SiO2");
    vproto.push_back("ABC_hP3_183_a_a_a;3;3;183;5;5;hP3;a,c/a,z1,z2,z3;-;AuCN;AuCN");
    vproto.push_back("AB_hP6_183_c_ab;2;6;183;4;5;hP6;a,c/a,z1,z2,z3;-;CrFe3NiSn5;CrFe3NiSn5");
    vproto.push_back("AB4C_hP72_184_d_4d_d;3;72;184;5;20;hP72;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Al[PO4];Al");
    vproto.push_back("A3BC_hP30_185_cd_c_ab;3;30;185;5;11;hP30;a,c/a,z1,z2,x3,z3,x4,z4,x5,y5,z5;-;KNiCl3;KNiCl3");
    vproto.push_back("A3B_hP24_185_ab2c_c;2;24;185;4;10;hP24;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Cu3P;Cu3P");
    vproto.push_back("A3B_hP8_185_c_a;2;8;185;4;5;hP8;a,c/a,z1,x2,z2;-;beta-RuCl3;beta-RuCl3");
    vproto.push_back("AB3_hP24_185_c_ab2c;2;24;185;4;10;hP24;a,c/a,z1,z2,x3,z3,x4,z4,x5,z5;-;Na3As;Na3As");
    vproto.push_back("A3B7_hP20_186_c_b2c;2;20;186;4;9;hP20;a,c/a,z1,x2,z2,x3,z3,x4,z4;D10_{2};Fe3Th7;Fe3Th7");
    vproto.push_back("AB3_hP4_187_e_fh;2;4;187;4;3;hP4;a,c/a,z3;-;Re3N;Re3N");
    vproto.push_back("A3BC_hP10_188_k_c_a;3;10;188;5;4;hP10;a,c/a,x3,y3;-;LiScI3;LiScI3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("A3BC_hP10_188_k_a_e;3;10;188;5;4;hP10;a,c/a,x3,y3;-;LiScI3;LiScI3");
    vproto.push_back("AB9C4_hP28_188_e_kl_ak;3;28;188;5;9;hP28;a,c/a,x3,y3,x4,y4,x5,y5,z5;S3_{2};BaSi4O9;BaSi4O9"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A8BC3D6_hP18_189_bfh_a_g_i;4;18;189;6;7;hP18;a,c/a,x3,x4,z5,x6,z6;E9_{b};pi-FeMg3Al8Si6;pi-FeMg3Al8Si6");
    vproto.push_back("A9BC3D5_hP18_189_fi_a_g_bh;4;18;189;6;7;hP18;a,c/a,x3,x4,z5,x6,z6;-;pi-FeMg3Al9Si5;pi-FeMg3Al9Si5");
    vproto.push_back("A2B_hP18_190_gh_bf;2;18;190;4;6;hP18;a,c/a,z2,x3,x4,y4;-;Li2Sb;Li2Sb");
    vproto.push_back("A5B3_hP16_190_bdh_g;2;16;190;4;5;hP16;a,c/a,x3,x4,y4;-;alpha-Sm3Ge5;alpha-Sm3Ge5");
    vproto.push_back("AB_hP24_190_i_afh;2;24;190;4;8;hP24;a,c/a,z2,x3,y3,x4,y4,z4;-;FeS;Troilite");
    vproto.push_back("A2B3C18D6_hP58_192_c_f_lm_l;4;58;192;6;9;hP58;a,c/a,x3,y3,x4,y4,x5,y5,z5;G3_{1};Be3Al2Si6O18;Beryl");
    vproto.push_back("AB2_hP72_192_m_j2kl;2;72;192;4;10;hP72;a,c/a,x1,x2,x3,x4,y4,x5,y5,z5;-;AlPO4;AlPO4");
    vproto.push_back("A5B3_hP16_193_dg_g;2;16;193;4;4;hP16;a,c/a,x2,x3;D8_{8};Mn5Si3;Mavlyanovite"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A3B_hP16_194_gh_ac;2;16;194;4;3;hP16;a,c/a,x4;D0_{24};Ni3Ti;Ni3Ti");
    vproto.push_back("A5B2_hP28_194_ahk_ch;2;28;194;4;6;hP28;a,c/a,x3,x4,x5,z5;D8_{11};Co2Al5;Co2Al5");
    vproto.push_back("A9B3C_hP26_194_hk_h_a;3;26;194;5;6;hP26;a,c/a,x2,x3,x4,z4;E9_{c};Al9Mn3Si;Al9Mn3Si");
    vproto.push_back("A12BC4_cP34_195_2j_ab_2e;3;34;195;5;9;cP34;a,x3,x4,x5,y5,z5,x6,y6,z6;-;PrRu4P12;PrRu4P12");
    vproto.push_back("A12B2C_cF60_196_h_bc_a;3;15;196;5;4;cF60;a,x4,y4,z4;-;Cu2Fe[CN]6;Cu2Fe6");
    //DX 20180925 - wrong space group, should be 210: vproto.push_back("A12B36CD12_cF488_196_2h_6h_ac_fgh;4;122;196;6;30;cF488;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;MgB12H12[H2O]12;MgB12H1212");
    vproto.push_back("ABC3_cP20_198_a_a_b;3;20;198;5;6;cP20;a,x1,x2,x3,y3,z3;G3;NaClO3;Sodium Chlorate");
    vproto.push_back("A2B11_cP39_200_f_aghij;2;39;200;4;7;cP39;a,x2,x3,x4,x5,y6,z6;D8_{c};Mg2Zn11;Mg2Zn11");
    vproto.push_back("AB3C_cP60_201_be_fh_g;3;60;201;5;7;cP60;a,x2,x3,x4,x5,y5,z5;-;KSbO3;KSbO3");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB3C_cP60_201_ce_fh_g;3;60;201;5;7;cP60;a,x2,x3,x4,x5,y5,z5;-;KSbO3;KSbO3");
    vproto.push_back("A6B6C_cF104_202_h_h_c;3;26;202;5;5;cF104;a,y2,z2,y3,z3;-;KB6H6;KB6H6");
    vproto.push_back("A_cF240_202_h2i;1;60;202;3;9;cF240;a,y1,z1,x2,y2,z2,x3,y3,z3;-;C;FCC C60 Buckminsterfullerine");
    vproto.push_back("A2BCD3E6_cF208_203_e_c_d_f_g;5;52;203;7;6;cF208;a,x3,x4,x5,y5,z5;-;Na3Co(CO3)2Cl;Pyrochlore");
    vproto.push_back("A4B2C6D16E_cF232_203_e_d_f_eg_a;5;58;203;7;7;cF232;a,x3,x4,x5,x6,y6,z6;-;Na6Mg2(SO4)(CO3)4;Tychite");
    vproto.push_back("AB3C16_cF160_203_a_bc_eg;3;40;203;5;5;cF160;a,x4,x5,y5,z5;-;Rb3AsSe16;Rb3AsSe16");//DX 20180925 - shifted Wyckoff positions, orig: vproto.push_back("AB3C16_cF160_203_b_ad_eg;3;40;203;5;5;cF160;a,x4,x5,y5,z5;-;Rb3AsSe16;Rb3AsSe16");
    vproto.push_back("A2B3C6_cP264_205_2d_ab2c2d_6d;3;264;205;5;33;cP264;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Ca3Al2O6;Ca3Al2O6");
    vproto.push_back("A_cP240_205_10d;1;240;205;3;31;cP240;a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;C;Simple Cubic C60 Buckminsterfullerine");
    vproto.push_back("AB3C2_cI96_206_c_e_ad;3;48;206;5;6;cI96;a,x2,x3,x4,y4,z4;E9_{d};AlLi3N2;AlLi3N2");
    vproto.push_back("A17B15_cP64_207_acfk_eij;2;64;207;4;8;cP64;a,x3,x4,y5,y6,x7,y7,z7;-;Pd17Se15;Pd17Se15");
    vproto.push_back("A3B_cP16_208_j_b;2;16;208;4;2;cP16;a,x2;-;PH3;PH3");
    vproto.push_back("A6B2CD6E_cP64_208_m_ad_b_m_c;5;64;208;7;7;cP64;a,x5,y5,z5,x6,y6,z6;-;Cs2ZnFe[CN]6;Cs2ZnFe6");
    vproto.push_back("A24BC_cF104_209_j_a_b;3;26;209;5;4;cF104;a,x3,y3,z3;-;F6KP;F6KP");
    vproto.push_back("A12B36CD12_cF488_210_h_3h_a_fg;4;122;210;6;15;cF488;a,x2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;MgB12H12[H2O]12;MgB12H1212"); //DX 20180925 - moved this structure from SG196 to SG210
    vproto.push_back("A12B6C_cF608_210_4h_2h_e;3;152;210;5;20;cF608;a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Te[OH]6;Te6");
    vproto.push_back("A2B_cI72_211_hi_i;2;36;211;4;4;cI72;a,y1,y2,y3;-;SiO2;SiO2");
    vproto.push_back("A2B_cP12_212_c_a;2;12;212;4;2;cP12;a,x2;-;SrSi2;SrSi2");
    vproto.push_back("A3B3C_cI56_214_g_h_a;3;28;214;5;3;cI56;a,y2,y3;-;Ca3PI3;Ca3PI3");
    vproto.push_back("A3BC2_cI48_214_f_a_e;3;24;214;5;3;cI48;a,x2,x3;-;Ag3AuTe2;Petzite");
    vproto.push_back("A4B9_cP52_215_ei_3efgi;2;52;215;4;11;cP52;a,x1,x2,x3,x4,x5,x6,x7,z7,x8,z8;D8_{3};gamma-Cu9Al4;gamma-brass");
    vproto.push_back("ABCD_cF16_216_c_d_b_a;4;4;216;6;1;cF16;a;-;LiMgAuSn;Quaternary Heusler"); //DX 20180925 - fixed typo in "Quaternary"
    vproto.push_back("A3B4C_cP16_218_c_e_a;3;16;218;5;2;cP16;a,x3;H2_{1};Ag3[PO4];Ag3"); //DX 20180925 - added Strukturbericht
    vproto.push_back("A7BC3D13_cF192_219_de_b_c_ah;4;48;219;6;5;cF192;a,x5,x6,y6,z6;-;Mg3B7ClO13;Boracite");
    vproto.push_back("A15B4_cI76_220_ae_c;2;38;220;4;5;cI76;a,x2,x3,y3,z3;D8_{6};Cu15Si4;Cu15Si4");
    vproto.push_back("A4B3_cI28_220_c_a;2;14;220;4;2;cI28;a,x2;D7_{3};Th3P4;Th3P4");
    vproto.push_back("A2B3C6_cP33_221_cd_ag_fh;3;33;221;5;4;cP33;a,x4,x5,x6;E9_{1};Ca3Al2O6;Ca3Al2O6");
    vproto.push_back("A5B3C16_cP96_222_ce_d_fi;3;96;222;5;6;cP96;a,x3,x4,x5,y5,z5;-;Ce5Mo3O16;Ce5Mo3O16");
    vproto.push_back("A23B6_cF116_225_bd2f_e;2;29;225;4;4;cF116;a,x3,x4,x5;D8_{a};Th6Mn23;Th6Mn23");
    vproto.push_back("A6B2C_cF36_225_e_c_a;3;9;225;5;2;cF36;a,x3;J1_{1};K2PtCl6;K2PtCl6");
    vproto.push_back("AB13_cF112_226_a_bi;2;28;226;4;3;cF112;a,y3,z3;D2_{3};NaZn13;NaZn13");
    vproto.push_back("A2B2C7_cF88_227_c_d_af;3;22;227;5;2;cF88;a,x4;-;Eu2Ir2O7;Pyrochlore Iridate");
    vproto.push_back("A3B4_cF56_227_ad_e;2;14;227;4;2;cF56;a,x3;D7_{2}/-;Co3O4/Co3O4;Spinel/Co3O4 (ICSD #63164)"); //DX20210106 - added metal-oxide
    vproto.push_back("A5BCD6_cF416_228_eg_c_b_h;4;104;228;6;6;cF416;a,x3,y4,x5,y5,z5;-;CuCrCl5[NH3]6;CuCrCl56");
    vproto.push_back("A6B_cF224_228_h_c;2;56;228;4;4;cF224;a,x2,y2,z2;-;TeO6H6;TeO6H6");
    vproto.push_back("A3B10_cI52_229_e_fh;2;26;229;4;4;cI52;a,x1,x2,y3;D8_{1};gamma-Fe3Zn10;gamma-brass");
    vproto.push_back("A4B_cI10_229_c_a;2;5;229;4;1;cI10;a;-;beta-Hg4Pt;beta-Hg4Pt");
    vproto.push_back("A7B3_cI40_229_df_e;2;20;229;4;3;cI40;a,x2,x3;D8_{f};Ir3Ge7;Ir3Ge7");
    vproto.push_back("A2B3C12D3_cI160_230_a_c_h_d;4;80;230;6;4;cI160;a,x4,y4,z4;S1_{4};Co3Al2Si3O12;Garnet");
    //DX 20181130 - add Ohad's SQS structures - START
    // -------------------------------------------------------------------------
    // SQS (from O. Levy)
    // -------------------------------------------------------------------------
    vproto.push_back("AB_aP16_2_4i_4i;2;16;2;4;30;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-/-;TaTi/Mo1O2;TaTi/Mo1O2 (ICSD #36263"); //DX20210105 - added metal-oxide
    vproto.push_back("A5B11_mP16_6_2abc_2a3b3c;2;16;6;4;32;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Ta5Ti11;Ta5Ti11");
    vproto.push_back("AB3_mC32_8_4a_12a;2;16;8;4;36;mC32;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,z15,x16,z16;-;TaTi3;TaTi3");
    vproto.push_back("AB3_mC32_8_4a_4a4b;2;16;8;4;32;mC32;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;TaTi3;TaTi3");
    vproto.push_back("A3B13_oC32_38_ac_a2bcdef;2;16;38;4;18;oC32;a,b/a,c/a,z1,z2,z3,z4,x5,z5,x6,z6,y7,z7,y8,z8,x9,y9,z9;-;Ta3Ti13;Ta3Ti13");
    vproto.push_back("A3B5_oC32_38_abce_abcdf;2;16;38;4;18;oC32;a,b/a,c/a,z1,z2,z3,z4,x5,z5,x6,z6,y7,z7,y8,z8,x9,y9,z9;-;Ta3Ti5;Ta3Ti5");
    vproto.push_back("AB7_hR16_166_c_c2h;2;16;166;4;8;hR16;a,c/a,x1,x2,x3,z3,x4,z4;-;TaTi7;TaTi7");
    //DX 20181130 - add Ohad's SQS structures - END
    //DX 20181211 - add Corey's kesterite structure - START
    vproto.push_back("A2BCD4_tI16_82_ac_b_d_g;4;8;82;6;5;tI16;a,c/a,x5,y5,z5;-;Cu2(Zn,Fe)SnS4;Kesterite");
    //DX 20181211 - add Corey's kesterite structure - END
    // -------------------------------------------------------------------------
    // misc prototypes (from Y. Lederer) //DX 20200203
    // -------------------------------------------------------------------------
    vproto.push_back("AB3_mC8_12_a_di;2;4;12;4;6;mC8;a,b/a,c/a,beta,x3,z3;-/-/-;-/Li1N3/Li1N3;-/Li1N3 (ICSD #34675)/Li1N3 (ICSD #181559)"); //DX20201019 - added two metal nitrides
    vproto.push_back("AB_mC8_12_i_i;2;4;12;4;8;mC8;a,b/a,c/a,beta,x1,z1,x2,z2;-;-;-");
    vproto.push_back("AB3C4_mC16_12_a_di_2i;3;8;12;5;10;mC16;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5;-;-;-");
    vproto.push_back("ABC2_mC16_12_i_i_adi;3;8;12;5;10;mC16;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5;-;-;-");
    vproto.push_back("AB3_oP4_47_a_ct;2;4;47;4;4;oP4;a,b/a,c/a,z3;-;-;-");
    vproto.push_back("A3B_oP4_47_cr_a;2;4;47;4;4;oP4;a,b/a,c/a,z3;-;-;-");
    vproto.push_back("AB3C4_oP8_47_a_ct_egs;3;8;47;5;5;oP8;a,b/a,c/a,z5,z6;-;-;-");
    vproto.push_back("A3BC4_oP8_47_eq_g_bdt;3;8;47;5;5;oP8;a,b/a,c/a,z5,z6;-;-;-");
    vproto.push_back("AB_oP4_51_e_e;2;4;51;4;5;oP4;a,b/a,c/a,z1,z2;-;-;-");
    vproto.push_back("ABC2_oP8_51_e_e_2f;3;8;51;5;7;oP8;a,b/a,c/a,z1,z2,z3,z4;-;-;-");
    vproto.push_back("AB_oP4_59_a_a;2;4;59;4;5;oP4;a,b/a,c/a,z1,z2;-;-;-");
    vproto.push_back("ABC2_oP8_59_a_a_2b;3;8;59;5;7;oP8;a,b/a,c/a,z1,z2,z3,z4;-/-/-;-/Mn1Na1O2/Li1Mn1O2;-/Mn1Na1O2 (ICSD #16271)/Li1Mn1O2 (ICSD #84642)"); //DX20210106 - added two metal-oxides
    vproto.push_back("ABC2_oC16_63_c_c_g;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,x3,y3;-;-;-");
    vproto.push_back("AB2C3_oC12_65_a_i_cj;3;6;65;5;5;oC12;a,b/a,c/a,y3,y4;-;-;-");
    vproto.push_back("A3BC4_oC16_65_ai_b_q;3;8;65;5;6;oC16;a,b/a,c/a,y3,x4,y4;-;-;-");
    vproto.push_back("A3BC4_oC16_65_bj_a_eh;3;8;65;5;5;oC16;a,b/a,c/a,x4,y5;-;-;-");
    vproto.push_back("AB3C4_oC16_65_a_bf_hi;3;8;65;5;5;oC16;a,b/a,c/a,x4,y5;-;-;-");
    vproto.push_back("AB2_oC6_65_a_i;2;3;65;4;4;oC6;a,b/a,c/a,y2;-;-2;-2");
    vproto.push_back("A3B_oC8_65_ai_b;2;4;65;4;4;oC8;a,b/a,c/a,y3;-;-;-");
    vproto.push_back("A3B_oC8_65_bj_a;2;4;65;4;4;oC8;a,b/a,c/a,y3;-;-;-");
    vproto.push_back("AB_oC8_65_i_i;2;4;65;4;5;oC8;a,b/a,c/a,y1,y2;-;-;-");
    vproto.push_back("ABC2_oC16_65_i_i_fh;3;8;65;5;6;oC16;a,b/a,c/a,x2,y3,y4;-;-;-");
    vproto.push_back("AB2C3_oI12_71_a_e_df;3;6;71;5;5;oI12;a,b/a,c/a,x3,x4;-;-;-");
    vproto.push_back("AB_tP2_123_a_b;2;2;123;4;2;tP2;a,c/a;-;-;-");
    vproto.push_back("AB_tP2_123_a_c;2;2;123;4;2;tP2;a,c/a;-;-;-");
    vproto.push_back("A2B_tP3_123_g_a;2;3;123;4;3;tP3;a,c/a,z2;-/-;-/O2Zr1;-/O2Zr1 (ICSD #92091)"); //DX20210106 - added metal-oxide
    vproto.push_back("A3B_tP4_123_abc_d;2;4;123;4;2;tP4;a,c/a;-;-;-");
    vproto.push_back("A3B_tP4_123_ag_b;2;4;123;4;3;tP4;a,c/a,z3;-;-;-");
    vproto.push_back("A3B_tP4_123_cf_a;2;4;123;4;2;tP4;a,c/a;-;-;-");
    vproto.push_back("AB3_tP4_123_a_bh;2;4;123;4;3;tP4;a,c/a,z3;-;-;-");
    vproto.push_back("ABC2_tP4_123_a_b_h;3;4;123;5;3;tP4;a,c/a,z3;-;-;-");
    vproto.push_back("ABC2_tP4_123_a_c_e;3;4;123;5;2;tP4;a,c/a;-/-;-/Fe1Ni1Pt2;-/Fe1Ni1Pt2 (ICSD #42564)"); //DX20201103 - added metal
    vproto.push_back("ABC2_tP4_123_a_d_bc;3;4;123;5;2;tP4;a,c/a;-;-;-");
    vproto.push_back("AB_tP4_123_g_g;2;4;123;4;4;tP4;a,c/a,z1,z2;-;-;-");
    vproto.push_back("A2BC3_tP6_123_g_b_ch;3;6;123;5;4;tP6;a,c/a,z3,z4;-;-;-");
    vproto.push_back("A3BC4_tP8_123_abc_d_i;3;8;123;5;3;tP8;a,c/a,z5;-;-;-");
    vproto.push_back("A3BC4_tP8_123_ag_b_2h;3;8;123;5;5;tP8;a,c/a,z3,z4,z5;-;-;-");
    vproto.push_back("A3BC4_tP8_123_cf_a_k;3;8;123;5;3;tP8;a,c/a,x4;-;-;-");
    vproto.push_back("AB3C4_tP8_123_a_bh_cdg;3;8;123;5;4;tP8;a,c/a,z5,z6;-;-;-");
    vproto.push_back("ABC2_tP8_123_h_h_abg;3;8;123;5;5;tP8;a,c/a,z3,z4,z5;-;-;-");
    vproto.push_back("ABC2_tP8_129_c_c_2c;3;8;129;5;6;tP8;a,c/a,z1,z2,z3,z4;-;-;-");
    vproto.push_back("A3B_tI8_139_ae_b;2;4;139;4;3;tI8;a,c/a,z3;-;-;-");
    vproto.push_back("AB3_tI8_139_a_bd;2;4;139;4;2;tI8;a,c/a;-;-;-");
    vproto.push_back("AB2C3_tI12_139_a_e_be;3;6;139;5;4;tI12;a,c/a,z3,z4;-;-;-");
    vproto.push_back("A3BC4_tI16_139_ae_b_g;3;8;139;5;4;tI16;a,c/a,z3,z4;-;-;-");
    vproto.push_back("AB3C4_tI16_139_a_bd_ce;3;8;139;5;3;tI16;a,c/a,z5;-;-;-");
    vproto.push_back("ABC2_tI16_139_e_e_cd;3;8;139;5;4;tI16;a,c/a,z3,z4;-;-;-");
    vproto.push_back("ABC2_tI16_141_a_b_e;3;8;141;5;3;tI16;a,c/a,z3;-/-;-/Fe2Li2O4;-/Fe2Li2O4 (ICSD #31149)"); //DX20210106 - added metal-oxide
    vproto.push_back("A2B_hP3_164_d_a;2;3;164;4;3;hP3;a,c/a,z2;-/-/-/-;-/-/N2Zr1/O2Pt1;-/-/N2Zr1 (ICSD #262746)/O2Pt1 (ICSD #76431)"); //DX20201103 -- added metal-nitride //DX20210106 - added metal-oxide
    vproto.push_back("A2BC3_hP6_164_d_a_bd;3;6;164;5;4;hP6;a,c/a,z3,z4;-/-;-/Bi2La1Li3;-/Bi2La1Li3 (ICSD #616769)"); //DX20201028 - added metal prototype
    vproto.push_back("A2BC3_hP6_164_d_b_ad;3;6;164;5;4;hP6;a,c/a,z3,z4;-;-;-");
    vproto.push_back("A3B_hR4_166_bc_a;2;4;166;4;3;hR4;a,c/a,x3;-/-;-/N3Na1;-/N3Na1 (ICSD #1144)"); //DX20201019 - added metal-nitride
    vproto.push_back("AB3_hR4_166_a_bc;2;4;166;4;3;hR4;a,c/a,x3;-;-;-");
    vproto.push_back("AB_hR4_166_c_c;2;4;166;4;4;hR4;a,c/a,x1,x2;-/-/-/-;-/-/B2Li2/C1Ru1;-/-/B2Li2 (ICSD #1)/C1Ru1 (ICSD #188285)"); //DX20210104 - added metal-boride //DX20210120 - added metal-carbide
    vproto.push_back("A3BC4_hR8_166_bc_a_2c;3;8;166;5;5;hR8;a,c/a,x3,x4,x5;-;-;-");
    vproto.push_back("AB3C4_hR8_166_a_bc_2c;3;8;166;5;5;hR8;a,c/a,x3,x4,x5;-;-;-");
    vproto.push_back("ABC2_hR8_166_c_c_abc;3;8;166;5;5;hR8;a,c/a,x3,x4,x5;-;-;-");
    vproto.push_back("AB3C4_cP8_221_a_c_bd;3;8;221;5;1;cP8;a;-;-;-");
    // -------------------------------------------------------------------------
    // oxide prototypes (from R. Friedrich) //DX20200624
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("AB_mC20_12_b2i_c2i;2;10;12;4;12;mC20;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6;-;OTi;OTi (ICSD #56694)");
    vproto.push_back("A2B3_mC20_12_2i_3i;2;10;12;4;14;mC20;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-/-;Ga2O3 (beta-Ga2O3)/Nb2O3;Ga2O3 (ICSD #83645) and beta-Ga2O3 (part 3)/Nb2O3 (ICSD #263119)"); //DX20210105 - added metal-oxide prototype //DX20210427 - added equivalent part 3 prototype
    vproto.push_back("A5B3_mC32_12_5i_3i;2;16;12;4;20;mC32;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-/-;O5Ti3/O5Ti3;O5Ti3 (ICSD #75193)/O5Ti3 (ICSD #26492)"); //DX20210105 - added metal-oxide prototype
    vproto.push_back("A3B_mP16_14_3e_e;2;16;14;4;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-/-/-;O3W/O3Rb1/O3Rb1/O3Rb1/O3W1;O3W (ICSD #84848)/O3Rb1 (ICSD #6094)/O3Rb1 (ICSD #47164)/O3Rb1 (ICSD #180568)/O3W1 (ICSD #647640)"); //DX20210105 - added four metal-oxide prototypes
    vproto.push_back("A2B3_mP20_14_2e_3e;2;20;14;4;19;mP20;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-/D5_{f};Bi2O3/As2S3;Bi2O3 (ICSD #168806)/Orpiment (As2S3, D5f)"); //DX20210427 - added Orpiment (As2S3, D5f) from part 3
    vproto.push_back("AB_mC8_15_a_e;2;4;15;4;5;mC8;a,b/a,c/a,beta,y2;-/-;CuO/Cu1O1;CuO (ICSD #92368)/Cu1O1 (ICSD #67850)"); //DX20210106 - added metal-oxide prototype
    vproto.push_back("A4B_mC20_15_2f_e;2;10;15;4;11;mC20;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3;-;O4Os;O4Os (ICSD #23803)");
    vproto.push_back("A2B5_oP28_19_2a_5a;2;28;19;4;24;oP28;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;As2O5;As2O5 (ICSD #654040)");
    vproto.push_back("A2B_oP24_33_4a_2a;2;24;33;4;21;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;O2Sb (O4Sb2);O2Sb (ICSD #4109) and alpha-Sb2O4 (cervantite, part 3)"); //DX20210428 - added equivalent part 3 prototype (alpha-Sb2O4, cervantite, http://aflow.org/prototype-encyclopedia/A2B_oP24_33_4a_2a.html)
    vproto.push_back("AB3_oC16_40_b_a2b;2;8;40;4;10;oC16;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-/-;CrO3/CrO3;CrO3 (ICSD #16031)/Orthorhombic CrO3"); //DX20200805 - created with bad geometry file; fixed //DX20210427 - added Orthorhombic CrO3 from part 3
    vproto.push_back("A5B2_oP14_59_a2e_e;2;14;59;4;10;oP14;a,b/a,c/a,z1,y2,z2,y3,z3,y4,z4;-/-;O5V2 (O5V2)/O5V2;O5V2 (ICSD #60767) and O5V2 (shcherbinaite (revised), part 3)/O5V2 (ICSD #43132)"); //DX20210106 - added metal-oxide prototype //DX20210428 - added equivalent part 3 prototype (O5V2, shcherbinaite (revised), http://aflow.org/prototype-encyclopedia/A5B2_oP14_59_a2f_f.html)
    vproto.push_back("AB_oC16_64_e_f;2;8;64;4;6;oC16;a,b/a,c/a,y1,y2,z2;-;KO;KO (ICSD #180559)");
    vproto.push_back("A4B3_tP28_135_gh_dh;2;28;135;4;7;tP28;a,c/a,x2,x3,y3,x4,y4;-;O4Pb3;O4Pb3:vs (ICSD #29094)");
    vproto.push_back("A2B3_hP15_144_2a_3a;2;15;144;4;17;hP15;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;B2O3;B2O3 (ICSD #24649)");
    vproto.push_back("A2B_hR3_166_c_a;2;3;166;4;3;hR3;a,c/a,x2;-/-/-/-;Cs2O/Ca2N1/Ba2N1/Li2O1;Cs2O (ICSD #27919)/Ca2N1 (ICSD #22231)/Ba2N1 (ICSD #409851)/Li2O1 (ICSD #108886)"); //DX20201019 - added two metal-nitrides //DX20210106 - added metal-oxide
    vproto.push_back("AB2_hR6_166_c_2c;2;6;166;4;5;hR6;a,c/a,x1,x2,x3;-/C12;OTl2/CaSi2;OTl2 (ICSD #16220)/CaSi2 (C12)"); //DX20210427 - added CaSi2 (C12) from part 3
    vproto.push_back("AB_hP12_189_fg_eh;2;12;189;4;6;hP12;a,c/a,z1,x2,x3,z4;-;NaO;NaO (ICSD #25526)");
    vproto.push_back("AB_hP8_194_ac_f;2;8;194;4;3;hP8;a,c/a,z3;-/-;LiO/N1Nb1;LiO (ICSD #152183)/N1Nb1 (ICSD #76008)"); //DX20201103 - added metal-nitride
    vproto.push_back("A2B3_cI80_199_a2b_2c;2;40;199;4;10;cI80;a,x1,x2,x3,x4,y4,z4,x5,y5,z5;-;In2O3;In2O3 (ICSD #33649)");
    vproto.push_back("A3B2_cF80_227_f_e;2;20;227;4;3;cF80;a,x1,x2;-/D6_{1};O3Sb2/O3Sb2;O3Sb2 (ICSD #31102)/Senarmontite (Sb2O3, D61)"); //DX20210427 - added Senarmontite (Sb2O3, D61) from part 3
    // ternaries
    vproto.push_back("A4B4C_aP18_2_4i_4i_i;3;18;2;5;33;aP18;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Na4O4Si;Na4O4Si (ICSD #15500)");
    vproto.push_back("A2B7C2_aP22_2_2i_7i_2i;3;22;2;5;39;aP22;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-/-;Ca2O7V2/Mg2O7V2;Ca2O7V2 (ICSD #421266)/Mg2O7V2 (ICSD #2321)");
    vproto.push_back("AB3C_aP30_2_3i_9i_3i;3;30;2;5;51;aP30;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;CaO3Si (CaO3Si);CaO3Si (ICSD #201537) and CaO3Si (wollastonite, part 3)"); //DX20210428 - added equivalent part 3 prototype (CaO3Si, wollastonite, http://aflow.org/prototype-encyclopedia/AB3C_aP30_2_3i_9i_3i.html)
    vproto.push_back("A2B5C_aP32_2_4i_10i_2i;3;32;2;5;54;aP32;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;S0_{1};Al2O5Si (Al2O5Si);Al2O5Si (ICSD #85742) and Al2O5Si (kyanite, part 3)"); //DX20210427 - equivalent to part 3 prototype
    vproto.push_back("A2B4C_mP28_4_4a_8a_2a;3;28;4;5;46;mP28;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Al2O4Sr;Al2O4Sr (ICSD #160297)");
    vproto.push_back("A3BC_mC60_5_ab8c_ab2c_3c;3;30;5;5;47;mC60;a,b/a,c/a,beta,y1,y2,y3,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;-;O3SiSr;O3SiSr (ICSD #32542)");
    vproto.push_back("A3B5C_mC54_8_3a3b_9a3b_3a;3;27;8;5;52;mC54;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;Ca3O5Si;Ca3O5Si (ICSD #81100)");
    vproto.push_back("A2B7C3_mP24_11_2e_7e_3e;3;24;11;5;28;mP24;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12;-;Na2O7Ti3;Na2O7Ti3 (ICSD #15463)");
    vproto.push_back("AB6C2_mC18_12_a_3i_i;3;9;12;5;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-/-/-/-;CaO6V2/MgO6V2/Co1O6V2/PdSe6Ta2;CaO6V2 (ICSD #166516)/MgO6V2 (ICSD #10391)/Co1O6V2 (ICSD #188911)/Ta2PdSe6"); //DX20210106 - add metal-oxide //DX20210427 - added Ta2PdSe6 from part 3
    vproto.push_back("ABC4_mC48_12_gi_hi_2i3j;3;24;12;5;23;mC48;a,b/a,c/a,beta,y1,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-/-;MgMoO4/FeMoO4;MgMoO4 (ICSD #20418)/FeMoO4 (ICSD #43012)");
    vproto.push_back("AB4C_mP12_13_f_2g_e;3;12;13;5;12;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;H0_{6}/-/-/-;MgO4W (MgO4W)/NiO4W/Mg1O4W1/Cu1O4W1;MgO4W (ICSD #67903) and MgO4W (huanzalaite, part 3)/NiO4W (ICSD #16685)/Mg1O4W1 (ICSD #36310)/Cu1O4W1 (ICSD #182751)"); //DX20210106 - added two metal-oxides //DX20210428 - added equivalent part 3 prototype (MgO4W, huanzalaite, H0_{6}, http://aflow.org/prototype-encyclopedia/AB4C_mP12_13_f_2g_e.html)
    vproto.push_back("AB3C_mP40_14_2e_6e_2e;3;40;14;5;34;mP40;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;MgO3Si;MgO3Si (ICSD #30895)");
    vproto.push_back("A2B3C_mC24_15_2e_af_e;3;12;15;5;10;mC24;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5;-;Li2O3Zr;Li2O3Zr (ICSD #94893)");
    vproto.push_back("AB3C_mC40_15_2e_3f_f;3;20;18;5;18;mC40;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;NaO3V;NaO3V (ICSD #2103)");
    vproto.push_back("A2B3C_mC48_15_aef_3f_2e;3;24;15;5;19;mC48;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Li2O3Ti (Na2O3Pr);Li2O3Ti (ICSD #162215) and Na2O3Pr (part 3)"); //DX20210427 - added equivalent part 3 prototype
    vproto.push_back("A4BC7_mC48_15_2f_e_e3f;3;24;15;5;21;mC48;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Al4CaO7;Al4CaO7 (ICSD #14270)");
    vproto.push_back("AB3C_mC60_15_cf_e4f_ef;3;30;15;5;24;mC60;a,b/a,c/a,beta,y2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;CaO3Si;CaO3Si (ICSD #87694)");
    vproto.push_back("ABC2_oP16_33_a_a_2a;3;16;33;5;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-/-/-/-;FeNaO2/Ga1Na1O2/Ga1Li1O2 (GaLiO2)/Fe1Na1O2/Ga1Na1O2/Ag1Al1O2;FeNaO2 (ICSD #186309)/Ga1Na1O2 (ICSD #4416)/Ga1Li1O2 (ICSD #18152) and GaLiO2 (part 3)/Fe1Na1O2 (ICSD #27117)/Ga1Na1O2 (ICSD #36652)/Ag1Al1O2 (ICSD #160643)"); //DX20210106 - added five metal-oxides //DX20210428 - added equivalent part 3 prototype (GaLiO2, http://aflow.org/prototype-encyclopedia/ABC2_oP16_33_a_a_2a.html)
    vproto.push_back("A2B3C_oC24_36_b_ab_a;3;12;36;5;13;oC24;a,b/a,c/a,y1,z1,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-;Li2O3Si/Na2O3Si/Cs2O3Pr1;Li2O3Si (ICSD #28192)/Na2O3Si (ICSD #24664)/Cs2O3Pr1 (ICSD #1181)"); //DX20210106 - added metal-oxide
    vproto.push_back("A2B5C_oP32_58_eg_3gh_g;3;32;58;5;17;oP32;a,b/a,c/a,z1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;SO_{2};Al2O5Si (Al2O5Si);Al2O5Si (ICSD #30679) and Al2O5Si (andalusite, part 3)"); //DX20210427 - added equivalent part 3 prototype info
    vproto.push_back("AB2C4_oC28_63_c_ac_fg;3;14;63;5;9;oC28;a,b/a,c/a,y2,y3,y4,z4,x5,y5;H1_{8};CrNa2O4 (CrNa2O4);CrNa2O4 (ICSD #76001) and CrNa2O4 (H1_{8}, part 3)"); //DX20210428 - added equivalent part 3 prototype (CrNa2O4, H1_{8}, http://aflow.org/prototype-encyclopedia/AB2C4_oC28_63_c_bc_fg.html)
    vproto.push_back("AB5C2_oC32_63_c_c2f_f;3;16;63;5;11;oC32;a,b/a,c/a,y1,y2,y3,z3,y4,z4,y5,z5;E4_{1}/-/-/-;MgO5Ti2 Fe2O5Ti/Mg1O5V2/Li1O5V2/NiS5Ta2;MgO5Ti2 (ICSD #37232) and Fe2O5Ti (pseudobrookite, part 3)/Mg1O5V2 (ICSD #50979)/Li1O5V2 (ICSD #50980)/Ta2NiS5"); //DX20210106 - added two metal-oxides //DX20210427 - added Ta2NiS5 from part 3 //DX20210427 - added Fe2O5Ti (E4_{1}, pseudobrookite, part 3) info
    vproto.push_back("ABC4_tI24_88_a_b_f;3;12;88;5;5;tI24;a,c/a,x3,y3,z3;-/H0_{4}/-;BaMoO4/CaMoO4 (CaO4W)/CaMoO4;BaMoO4 (ICSD #50821)/CaMoO4 (ICSD #417513) and Scheelite (CaO4W, part 3)/CaMoO4 (ICSD #77334)"); //DX20210428 - added equivalent part 3 prototype (CaO4W, scheelite, H0_{4}, http://aflow.org/prototype-encyclopedia/AB4C_tI24_88_b_f_a.html)
    vproto.push_back("A4BC2_tP28_91_2d_b_ac;3;28;91;5;11;tP28;a,c/a,y1,y2,x3,x4,y4,z4,x5,y5,z5;-;O4TiZn2;O4TiZn2 (ICSD #109093)");
    vproto.push_back("A2BC4_hP56_173_2b2c_ac_b5c;3;56;173;5;30;hP56;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;Al2BaO4;Al2BaO4 (ICSD #246027)");
    vproto.push_back("AB2C4_cF56_227_b_c_e;3;14;227;5;2;cF56;a,x3;-;MoNa2O4;MoNa2O4 (ICSD #151970)");
    // -------------------------------------------------------------------------
    // nitrides prototypes (from R. Friedrich) //DX20200624
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A2B3_cI80_206_ad_e;2;40;206;4;5;cI80;a,x2,x3,y3,z3;-;N2Zn3;N2Zn3 (ICSD #84918)");
    // ternaries
    vproto.push_back("A3B3C_mC56_9_6a_6a_2a;3;28;9;5;46;mC56;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;N3Na3W1;N3Na3W1 (ICSD #75364)"); //DX20210604 - added nitride for R. Friedrich (N3Na3W1, ICSD #75364)
    vproto.push_back("ABC_oP12_62_c_c_c;3;12;62;5;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-;CaLiN/Li1Mg1N1/Ca1Mg1Sn1 (CuMnP)/Ni1Sn1Tb1/Co1Dy1Sn1/Bi1Ca1Li1/Ni1Sc1Sn1/Co1Dy1Sn1/La1Ni1Sn1/Al1Ce1Pt1/La1Ni1Sn1/Pt1Yb1Zn1/Er1Ga1Rh1/Ga1Pd1Sc1/Li1O1Rb1/C1Li1N1;CaLiN (ICSD #107304)/Li1Mg1N1 (ICSD #93259)/Ca1Mg1Sn1 (ICSD #42757) and CuMnP (part 3)/Ni1Sn1Tb1 (ICSD #54301)/Co1Dy1Sn1 (ICSD #54415)/Bi1Ca1Li1 (ICSD #58762)/Ni1Sc1Sn1 (ICSD #105338)/Co1Dy1Sn1 (ICSD #106459)/La1Ni1Sn1 (ICSD #108571)/Al1Ce1Pt1 (ICSD #150172)/La1Ni1Sn1 (ICSD #157921)/Pt1Yb1Zn1 (ICSD #159305)/Er1Ga1Rh1 (ICSD #630573)/Ga1Pd1Sc1 (ICSD #635080)/C1Li1N1 (ICSD #77321)"); //DX20201016 - added another metal-nitride system //DX20201028 - added 12 metal systems //DX20210106 - added metal-oxide //DX20210221 - added carbo-nitride //DX20210428 - added equivalent part 3 prototype (CuMnP, http://aflow.org/prototype-encyclopedia/ABC_oP12_62_c_c_c.html)
    vproto.push_back("A2B2C_tI10_139_e_e_a;3;5;139;5;4;tI10;a,c/a,z2,z3;-/-/-;Ca2N2Zn/N2Sr2Zn1/Ba2N2Zn1;Ca2N2Zn (ICSD #69049)/N2Sr2Zn1 (ICSD #80376)/Ba2N2Zn1 (ICSD #80377)"); //DX20201016 - added two metal-nitride prototypes
    vproto.push_back("A3BC2_hP6_191_f_a_d;3;6;191;5;2;hP6;a,c/a;-;Na3TaTi2;Na3TaTi2 (ICSD #186418)");

#if !(USE_HARDCODED_PROTOTYPES) //DX20210114 - these new prototypes are not hard-coded
    // -------------------------------------------------------------------------
    // metal-nitride prototypes (from DX) //DX20201016
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A2B_mP6_6_3ab_ab;2;6;6;4;16;mP6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;N2Re1;N2Re1 (ICSD #187448)");
    vproto.push_back("A2B_mC12_8_4a_2a;2;6;8;4;16;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;N2Re1;N2Re1 (ICSD #187447)");
    vproto.push_back("A3B_mC16_8_6a_b;2;8;8;4;19;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;N3Rb1;N3Rb1 (ICSD #155169)");
    vproto.push_back("A2B_mP6_11_2e_e;2;6;11;4;10;mP6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;N2Re1;N2Re1 (ICSD #187444)");
    vproto.push_back("AB6_mP14_11_e_6e;2;14;11;4;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Ba1N6;Ba1N6 (ICSD #14244)");
    vproto.push_back("A3B_mC8_12_di_a;2;4;12;4;6;mC8;a,b/a,c/a,beta,x3,z3;-/-;N3Na1/N3Na1;N3Na1 (ICSD #29370)/N3Na1 (ICSD #29376)");
    vproto.push_back("A2B_mC12_12_2i_i;2;6;12;4;10;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-/-/-/-/-;N2Re1/B2Mo1/C2Ca1/C2Ca1 (CaC2)/C2Ir1;N2Re1 (ICSD #187441)/B2Mo1 (ICSD #418398)/C2Ca1 (ICSD #54185)/C2Ca1 (ICSD #54188) and CaC2 (part 3)/C2Ir1 (ICSD #181488)"); //DX20210104 - added metal-boride //DX20210120 //DX20210428 - added equivalent part 3 prototype (CaC2-III, http://aflow.org/prototype-encyclopedia/A2B_mC12_12_2i_i.html)
    vproto.push_back("AB_mC16_12_2i_2i;2;8;12;4;12;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;N1Sr1;N1Sr1 (ICSD #411555)");
    vproto.push_back("A3B2_mC30_12_a4i_3i;2;15;12;4;18;mC30;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-/-;Ca3N2/Ca3N2;Ca3N2 (ICSD #162794)/Ca3N2 (ICSD #169726)");
    vproto.push_back("A2B_mP6_13_g_e;2;6;13;4;8;mP6;a,b/a,c/a,beta,y1,x2,y2,z2;-;N2Re1;N2Re1 (ICSD #187451)");
    vproto.push_back("AB2_mP12_14_e_2e;2;12;14;4;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ir1N2;Ir1N2 (ICSD #160623)");
    vproto.push_back("AB3_mP16_14_e_3e;2;16;14;4;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-;Hg2N6/B2H6/B2H6;Hg2N6 (ICSD #98661)/beta-B2H6/B2H6 (P21/c)"); //DX20210427 - added beta-B2H6 and B2H6 (P21/c) from part 3
    vproto.push_back("AB2_mC12_15_e_f;2;6;15;4;8;mC12;a,b/a,c/a,beta,y1,x2,y2,z2;C_{g};Ba1N2 (C2Th);Ba1N2 (ICSD #280681) and C2Th (C_{g}, part 3)"); //DX20210428 - added equivalent part 3 prototype (C2Th, C_{g}, http://aflow.org/prototype-encyclopedia/A2B_mC12_15_f_e.html)
    vproto.push_back("AB_oC16_36_2a_2a;2;8;36;4;11;oC16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-/-;N1Os1/O1Sn1;N1Os1 (ICSD #167514)/O1Sn1 (ICSD #424729)"); //DX20210105 - added metal-oxide prototype
    vproto.push_back("AB_oP8_53_h_h;2;8;53;4;7;oP8;a,b/a,c/a,y1,z1,y2,z2;-;N1Os1;N1Os1 (ICSD #167513)");
    vproto.push_back("A2B_oP12_57_e_c;2;12;57;4;7;oP12;a,b/a,c/a,x1,x2,y2,z2;-;N2Re1;N2Re1 (ICSD #187450)");
    vproto.push_back("A2B_oP12_60_d_c;2;12;60;4;7;oP12;a,b/a,c/a,y1,x2,y2,z2;-/-/-/-;Fe2N1 (Fe2N)/O2Ti1/O2Pb1 (O2Pb)/O2Re1;Fe2N1 (ICSD #150889) and zeta-Fe2N (part 3)/O2Ti1 (ICSD #15328)/O2Pb1 (ICSD #20362) and alpha-PbO2 (part 3)/O2Re1 (ICSD #24060)"); //DX20210106 - added four metal-oxide prototypes //DX20210428 - added equivalent part 3 prototype (zeta-Fe2N, http://aflow.org/prototype-encyclopedia/A2B_oP12_60_d_c.Fe2N.html) //DX20210428 - added equivalent part 3 prototype (alpha-PbO2, http://aflow.org/prototype-encyclopedia/A2B_oP12_60_d_c.html)
    vproto.push_back("A5B3_oC32_63_c2f_cf;2;16;63;4;11;oC32;a,b/a,c/a,y1,y2,y3,z3,y4,z4,y5,z5;-/-;N5Ta3/N5Ta3;N5Ta3 (ICSD #16253)/N5Ta3 (ICSD #66533)");
    vproto.push_back("A2B_oC12_65_gj_e;2;6;65;4;5;oC12;a,b/a,c/a,x2,y3;-;N2Re1;N2Re1 (ICSD #187452)");
    vproto.push_back("AB_oI4_71_b_c;2;4;71;4;3;oI4;a,b/a,c/a;-;Cr1N1;Cr1N1 (ICSD #53146)");
    vproto.push_back("A2B_oI12_71_2e_f;2;6;71;4;6;oI12;a,b/a,c/a,x1,x2,x3;-;N2Re1;N2Re1 (ICSD #187449)");
    vproto.push_back("AB_oI16_71_abe_n;2;8;71;4;6;oI16;a,b/a,c/a,x3,x4,y4;-;Li2N2;Li2N2 (ICSD #423831)");
    vproto.push_back("AB3_oI16_72_b_cj;2;8;72;4;5;oI16;a,b/a,c/a,x3,y3;-;Ag1N3;Ag1N3 (ICSD #27135)");
    vproto.push_back("A5B4_tI18_87_ah_h;2;9;87;4;6;tI18;a,c/a,x2,y2,x3,y3;-;N5Nb4;N5Nb4 (ICSD #26251)");
    vproto.push_back("AB3_tI32_88_d_cf;2;16;88;4;5;tI32;a,c/a,x3,y3,z3;-;Cu1N3 (CuN3);Cu1N3 (ICSD #30633) and Copper (I) Azide (CuN3, part 3)"); //DX20210428 - added equivalent part 3 prototype (Copper (I) Azide, CuN3, http://aflow.org/prototype-encyclopedia/AB3_tI32_88_d_cf.html)
    vproto.push_back("AB_tP2_123_d_a;2;2;123;4;2;tP2;a,c/a;-;N1Pr1;N1Pr1 (ICSD #168645)");
    vproto.push_back("A3B_tP4_123_ag_d;2;4;123;4;3;tP4;a,c/a,z3;-;N3Rb1;N3Rb1 (ICSD #16963)");
    vproto.push_back("A3B_tP8_123_egh_ab;2;8;123;4;4;tP8;a,c/a,z4,z5;-;Cu3N1;Cu3N1 (ICSD #180237)");
    vproto.push_back("A2B_tP6_127_g_b;2;6;127;4;3;tP6;a,c/a,x2;-;N2Re1;N2Re1 (ICSD #187442)");
    vproto.push_back("AB_tP4_131_e_c;2;4;131;4;2;tP4;a,c/a;-;N1Tc1;N1Tc1 (ICSD #187709)");
    vproto.push_back("AB_tI4_139_a_b;2;2;139;4;2;tI4;a,c/a;-/- (L2_{0});Mn1N1/Co1O1 (FeC_{x});Mn1N1 (ICSD #106932)/Co1O1 (ICSD #174027) and \"Martensite Type\" FeC_{x} (part 3)"); //DX20210106 - added metal-oxide //DX20210428 - added equivalent part 3 prototype ("Martensite Type" FeC_{x}, L2_{0}, http://aflow.org/prototype-encyclopedia/AB_tI4_139_b_a.html, part 3)
    vproto.push_back("A3B2_tI10_139_ae_e;2;5;139;4;4;tI10;a,c/a,z2,z3;-/-;Mn3N2/Mn3N2;Mn3N2 (ICSD #71638)/Mn3N2 (ICSD #84202)");
    vproto.push_back("A3B4_tI14_139_ad_ce;2;7;139;4;3;tI14;a,c/a,z4;-;N3Nb4;N3Nb4 (ICSD #76389)");
    vproto.push_back("A8B_tI18_139_deh_a;2;9;139;4;4;tI18;a,c/a,z3,x4;D2_{g};Fe16N2 (Fe8N);Fe16N2 (ICSD #41953) and Fe8N (part 3)"); //DX20210428 - added equivalent part 3 prototype (Fe8N, D2_{g}, http://aflow.org/prototype-encyclopedia/A8B_tI18_139_deh_a.html)
    vproto.push_back("AB3_tI16_140_a_dh;2;8;140;4;3;tI16;a,c/a,x3;-/-/-;K1N3/Cs1N3/Ag1N3;K1N3 (ICSD #24007)/Cs1N3 (ICSD #25008)/Ag1N3 (ICSD #183201)");
    vproto.push_back("AB3_tI16_140_a_ch;2;8;140;4;3;tI16;a,c/a,x3;-;Cu1N3;Cu1N3 (ICSD #187018)");
    vproto.push_back("AB3_tI16_140_c_ah;2;8;140;4;3;tI16;a,c/a,x3;-;K1N3;K1N3 (ICSD #1145)");
    vproto.push_back("AB2_tI12_141_a_e;2;6;141;4;3;tI12;a,c/a,z2;-;N1Ti2;N1Ti2 (ICSD #23403)");
    vproto.push_back("A3B_hR4_160_3a_a;2;4;160;4;6;hR4;a,c/a,x1,x2,x3,x4;-;N3Na1;N3Na1 (ICSD #644523)");
    vproto.push_back("A3B2_hR5_160_3a_2a;2;5;160;4;7;hR5;a,c/a,x1,x2,x3,x4,x5;-;Be3N2;Be3N2 (ICSD #185490)");
    vproto.push_back("AB_hP16_162_ek_ci;2;16;162;4;6;hP16;a,c/a,z2,x3,x4,z4;-;Mo1N1;Mo1N1 (ICSD #43559)");
    vproto.push_back("A3B2_hP5_164_ac_d;2;5;164;4;4;hP5;a,c/a,z2,z3;-;N3V2;N3V2 (ICSD #182700)");
    vproto.push_back("AB_hP16_164_ci_di;2;16;164;4;8;hP16;a,c/a,z1,z2,x3,z3,x4,z4;-;Mo1N1;Mo1N1 (ICSD #60168)");
    vproto.push_back("AB2_hR3_166_a_c;2;3;166;4;3;hR3;a,c/a,x2;-/-/-;N1Sr2/C1Y2/Ba1C2;N1Sr2 (ICSD #23530)/C1Y2 (ICSD #22283)/Ba1C2 (ICSD #186576)"); //DX20210120 - added two metal-carbides
    vproto.push_back("AB_hP16_186_ac_bc;2;16;186;4;8;hP16;a,c/a,z1,z2,x3,z3,x4,z4;-;Mo1N1;Mo1N1 (ICSD #76280)");
    vproto.push_back("AB2_hP3_187_a_h;2;3;187;4;3;hP3;a,c/a,z2;-;Hf1N2;Hf1N2 (ICSD #290427)");
    vproto.push_back("A2B_hP3_187_h_a;2;3;187;4;3;hP3;a,c/a,z2;-;N2W1;N2W1 (ICSD #290433)");
    vproto.push_back("AB_hP6_189_f_ad;2;6;189;4;3;hP6;a,c/a,x3;-;N1Ta1;N1Ta1 (ICSD #1396)");
    vproto.push_back("A2B_hP3_191_e_a;2;3;191;4;3;hP3;a,c/a,z2;-;N2Os1;N2Os1 (ICSD #260545)");
    vproto.push_back("A2B_hP6_191_af_d;2;6;191;4;2;hP6;a,c/a;-;N4Ta2;N4Ta2 (ICSD #182351)");
    vproto.push_back("AB_hP6_191_bc_ad;2;6;191;4;2;hP6;a,c/a;-;N1Y1;N1Y1 (ICSD #161079)");
    vproto.push_back("A3B_hP8_193_g_b;2;8;193;4;3;hP8;a,c/a,x2;-/-;Ba3N1/Cs3O1;Ba3N1 (ICSD #77730)/Cs3O1 (ICSD #15695)"); //DX20210106 - added metal-oxide
    vproto.push_back("AB_hP4_194_a_c;2;4;194;4;2;hP4;a,c/a;-/-/C_{k};Cd1N1/B1Pt1/LiZn2;Cd1N1 (ICSD #185566)/B1Pt1 (ICSD #24363)/LiZn2 (Ck)"); //DX20210104 - added metal-boride //DX20210427 - added LiZn2 (Ck)
    vproto.push_back("AB2_hP6_194_c_e;2;6;194;4;3;hP6;a,c/a,z2;-;Hf1N2;Hf1N2 (ICSD #290428)");
    vproto.push_back("A2B_hP6_194_e_c;2;6;194;4;3;hP6;a,c/a,z2;-/-;N2Ta1/C2Os1;N2Ta1 (ICSD #290431)/C2Os1 (ICSD #168280)"); //DX20210120 - added metal-carbide
    vproto.push_back("A3B2_hP10_194_bf_ac;2;10;194;4;3;hP10;a,c/a,z4;-/-;Be3N2/Ca3N2;Be3N2 (ICSD #25656)/Ca3N2 (ICSD #162797)");
    vproto.push_back("A3B2_hP10_194_cf_f;2;10;194;4;4;hP10;a,c/a,z2,z3;-/-;N3W2/B3Ru2;N3W2 (ICSD #186207)/B3Ru2 (ICSD #23715)"); //DX20210104 - added metal-boride
    vproto.push_back("AB_hP16_194_bh_ag;2;16;194;4;3;hP16;a,c/a,x4;-;Mo1N1;Mo1N1 (ICSD #106926)");
    vproto.push_back("AB3_cI32_204_c_g;2;16;204;4;3;cI32;a,y2,z2;-;Co1N3/Co1N3;Co1N3 (ICSD #162105)/Co1N3 (ICSD #162106)");
    vproto.push_back("A4B_cP5_221_ac_b;2;5;221;4;1;cP5;a;-;Mn4N1;Mn4N1 (ICSD #44369)");
    vproto.push_back("AB4_cP5_221_a_bc;2;5;221;4;1;cP5;a;L'1_{0};N1Ni4 (Fe4N);N1Ni4 (ICSD #76403) and Fe4N (part 3)"); //DX20210428 - added equivalent part 3 prototype (Fe4N, http://aflow.org/prototype-encyclopedia/A4B_cP5_221_bc_a.html)
    vproto.push_back("A4B3_cP7_221_ac_d;2;7;221;4;1;cP7;a;-;N4W3;N4W3 (ICSD #30370)");

    // -------------------------------------------------------------------------
    // ternaries
    vproto.push_back("AB2C5_aP8_1_a_2a_5a;3;8;1;5;30;aP8;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Li1N2Na5;Li1N2Na5 (ICSD #92316)");
    vproto.push_back("AB2C5_mC16_5_a_c_a2bc;3;8;5;5;14;mC16;a,b/a,c/a,beta,y1,y2,y3,y4,x5,y5,z5,x6,y6,z6;-;Li1N2Na5;Li1N2Na5 (ICSD #92315)");
    vproto.push_back("A3B2C3_mP8_6_2ab_ab_a2b;3;8;6;5;20;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;Li3N2Na3;Li3N2Na3 (ICSD #92312)");
    vproto.push_back("A5B2C_mP8_6_3a2b_ab_b;3;8;6;5;20;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;Li5N2Na1;Li5N2Na1 (ICSD #92314)");
    vproto.push_back("A2BC2_mC20_12_2i_i_2i;3;10;12;5;14;mC20;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Ca2Fe1N2;Ca2Fe1N2 (ICSD #72389)");
    vproto.push_back("AB2C12_mC30_12_a_i_6i;3;15;12;5;18;mC30;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-/-;Cd1K2N12/Cd1K2N12;Cd1K2N12 (ICSD #31297)/Cd1K2N12 (ICSD #659621)");
    vproto.push_back("ABC_mP12_14_e_e_e;3;12;14;5;13;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;-/E0_{7};Be1Li1N1/AsFeS;Be1Li1N1 (ICSD #402341)/Arsenopyrite (FeAsS, E07)"); //DX20210427 - added arsenopyrite from part 3
    vproto.push_back("A2BC_oC16_63_2c_c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-;Ca2In1N1;Ca2In1N1 (ICSD #96228)");
    vproto.push_back("ABC3_oC20_63_c_a_cf;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-/-;Al1N1Zr3/Ca1Ir1O3;Al1N1Zr3 (ICSD #29521)/Ca1Ir1O3 (ICSD #159026)"); //DX20210106 - added metal-oxide
    vproto.push_back("A3BC3_oC28_63_cg_c_cg;3;14;63;5;10;oC28;a,b/a,c/a,y1,y2,y3,x4,y4,x5,y5;-;Ca3Cr1N3;Ca3Cr1N3 (ICSD #40205)");
    vproto.push_back("A2BC_oC8_65_i_c_b;3;4;65;5;4;oC8;a,b/a,c/a,y3;-;Li4N2Na2;Li4N2Na2 (ICSD #92307)");
    vproto.push_back("A7B2C3_oC24_65_a2ij_i_cj;3;12;65;5;8;oC24;a,b/a,c/a,y3,y4,y5,y6,y7;-;Ca7N2Tl3;Ca7N2Tl3 (ICSD #98178)");
    vproto.push_back("A2BC3_oI12_71_e_a_de;3;6;71;5;5;oI12;a,b/a,c/a,x3,x4;-/-/-;Ce2Mn1N3/Ba2Cu1O3/Ca2Cu1O3;Ce2Mn1N3 (ICSD #50579)/Ba2Cu1O3 (ICSD #68217)/Ca2Cu1O3 (ICSD #93651"); //DX20210106 - added two metal-oxides
    vproto.push_back("ABC_tP6_99_ab_ab_ab;3;6;99;5;8;tP6;a,c/a,z1,z2,z3,z4,z5,z6;-;Cr1N1Nb1;Cr1N1Nb1 (ICSD #23779)");
    vproto.push_back("A3B2C3_tP8_115_abc_g_de;3;8;115;5;4;tP8;a,c/a,z5,z6;-;Li3N2Na3;Li3N2Na3 (ICSD #92311)");
    vproto.push_back("ABC_tP3_123_c_a_b;3;3;123;5;2;tP3;a,c/a;-;Fe1N1Ni1;Fe1N1Ni1 (ICSD #53505)");
    vproto.push_back("A5B2C_tP8_123_ai_bc_d;3;8;123;5;3;tP8;a,c/a,z5;-;Li5N2Na1;Li5N2Na1 (ICSD #92313)");
    vproto.push_back("ABC_tP6_129_c_c_c;3;6;129;5;5;tP6;a,c/a,z1,z2,z3;-/-;Ca1Ga1N1/Mo1N1Ta1;Ca1Ga1N1 (ICSD #2027)/Mo1N1Ta1 (ICSD #100437)");
    vproto.push_back("ABC2_tP8_129_c_c_bc;3;8;129;5;5;tP8;a,c/a,z2,z3,z4;-;Ba1Hf1N2;Ba1Hf1N2 (ICSD #50994)");
    vproto.push_back("ABC2_tP8_129_a_c_bc;3;8;129;5;4;tP8;a,c/a,z3,z4;-;Li2N2Na4;Li2N2Na4 (ICSD #92309)");
    vproto.push_back("ABC_tP6_131_e_c_b;3;6;131;5;2;tP6;a,c/a;-;Ca1N1Ni1;Ca1N1Ni1 (ICSD #69044)");
    vproto.push_back("ABC_tP6_131_b_c_e;3;6;131;5;2;tP6;a,c/a;-;Li1N1Sr1;Li1N1Sr1 (ICSD #87414)");
    vproto.push_back("A2BC_tI8_139_d_b_a;3;4;139;5;2;tI8;a,c/a;-;Li4N2Na2;Li4N2Na2 (ICSD #92305)");
    vproto.push_back("A3B6C2_tI22_139_ae_eg_e;3;11;139;5;6;tI22;a,c/a,z2,z3,z4,z5;-/-;La3N6V2/La3N6Nb2;La3N6V2 (ICSD #98477)/La3N6Nb2 (ICSD #411473)");
    vproto.push_back("A2B2C_tI20_140_h_h_a;3;10;140;5;4;tI20;a,c/a,x2,x3;-;Be4N4Sr2;Be4N4Sr2 (ICSD #413356)");
    vproto.push_back("A2BC2_tI20_140_h_a_h;3;10;140;5;4;tI20;a,c/a,x2,x3;-;Be4Ca2N4;Be4Ca2N4 (ICSD #413357)");
    vproto.push_back("AB2C2_tI20_140_a_h_h;3;10;140;5;4;tI20;a,c/a,x2,x3;-;Ba1Be2N2;Ba1Be2N2 (ICSD #415304)");
    vproto.push_back("AB2C_hR4_160_a_2a_a;3;4;160;5;6;hR4;a,c/a,x1,x2,x3,x4;-;Cr1N2W1;Cr1N2W1 (ICSD #84639)");
    vproto.push_back("A2B2C_hP5_164_d_d_a;3;5;164;5;4;hP5;a,c/a,z2,z3;-/-;Li2N2Zr1/Mg2N2Sr1 (and Ce2O2S);Li2N2Zr1 (ICSD #16231)/Mg2N2Sr1 (ICSD #410826) and Ce2O2S (part 3)");
    vproto.push_back("AB2C2_hP5_164_a_d_d;3;5;164;5;4;hP5;a,c/a,z2,z3;-;Ce1Li2N2;Ce1Li2N2 (ICSD #34003)");
    vproto.push_back("AB2C_hR4_166_a_c_b;3;4;166;5;3;hR4;a,c/a,x3;-/-/-/-;Cu1N2Ta1/Ba1Cu2Ga1/Hg1O2Sr1/C1N2Sr1;Cu1N2Ta1 (ICSD #71136)/Ba1Cu2Ga1 (ICSD #615828)/Hg1O2Sr1 (ICSD #165068)/C1N2Sr1 (ICSD #59860)"); //DX20201028 - added metal prototype //DX20210106 - added metal-oxide //DX20210221 - added carbo-nitride
    vproto.push_back("A3BC3_hP14_176_h_c_h;3;14;176;5;6;hP14;a,c/a,x2,y2,x3,y3;-;Ba3Fe1N3;Ba3Fe1N3 (ICSD #36502)");
    vproto.push_back("ABC_hP6_186_b_b_a;3;6;186;5;5;hP6;a,c/a,z1,z2,z3;-/-;N1Na1Sn1/Cu1Sn1Ti1;N1Na1Sn1 (ICSD #172471)/Cu1Sn1Ti1 (ICSD #54657)"); //DX20201028 - added two metal prototypes
    vproto.push_back("AB2C_hP8_186_a_ab_b;3;8;186;5;6;hP8;a,c/a,z1,z2,z3,z4;-;Mn1N2W1;Mn1N2W1 (ICSD #80029)");
    vproto.push_back("ABC_hP3_187_a_c_d;3;3;187;5;2;hP3;a,c/a;-;Li1N1Ni1;Li1N1Ni1 (ICSD #247028)");
    vproto.push_back("A5B3C3_hP11_189_dg_g_f;3;11;189;5;5;hP11;a,c/a,x2,x3,x4;-;Li5N3Ni3;Li5N3Ni3 (ICSD #411152)");
    vproto.push_back("A2BC_hP4_191_c_a_b;3;4;191;5;2;hP4;a,c/a;-;Li2N1Na1;Li2N1Na1 (ICSD #92308)");
    vproto.push_back("ABC2_hP4_191_b_a_c;3;4;191;5;2;hP4;a,c/a;-;Li1N1Na2;Li1N1Na2 (ICSD #92310)");
    vproto.push_back("ABC2_hP8_194_c_a_f;3;8;194;5;3;hP8;a,c/a,z3;-/-/-/-/-/-;Ba1Ce1N2/Ag1Fe1O2 (AlCuO2)/Cs1Nd1O2/Cu1Fe1O2/Al1C1W2/Al1C1Ta2;Ba1Ce1N2 (ICSD #74791)/Ag1Fe1O2 (ICSD #2786) and Hexagonal Delafossite (AlCuO2, part 3)/Cs1Nd1O2 (ICSD #27336)/Cu1Fe1O2 (ICSD #66546)/Al1C1W2 (ICSD #165101)/Al1C1Ta2 (ICSD #181247)"); //DX20210106 - added three metal-oxides //DX20210121 - added two metal-carbides //DX20210428 - added equivalent part 3 prototype (hexagonal delafossite, AlCuO2, http://aflow.org/prototype-encyclopedia/ABC2_hP8_194_a_c_f.html)
    vproto.push_back("AB2C_hP8_194_a_f_b;3;8;194;5;3;hP8;a,c/a,z3;-;Fe1N2W1;Fe1N2W1 (ICSD #75971)");
    vproto.push_back("ABC2_hP8_194_a_c_f;3;8;194;5;3;hP8;a,c/a,z3;-/-/-/-;Mg1Mo1N2/Na1Nb1O2/Al1Au1O2/C1Pb1Ti2;Mg1Mo1N2 (ICSD #185913)/Na1Nb1O2 (ICSD #29282)/Al1Au1O2 (ICSD #95663)/C1Pb1Ti2 (ICSD #42926)"); //DX20210106 - added three metal-oxides //DX20210121 - added metal-carbide
    vproto.push_back("A3BC_hP10_194_h_a_c;3;10;194;5;3;hP10;a,c/a,x3;-;Ba3N1Na1;Ba3N1Na1 (ICSD #67497)");
    vproto.push_back("AB3C4_hP16_194_c_ae_2f;3;16;194;5;5;hP16;a,c/a,z3,z4,z5;-/-/-;Al1N3V4/Al1C3Ta4/Al1C3Ta4;Al1N3V4 (ICSD #181351)/Al1C3Ta4 (ICSD #157843)/Al1C3Ta4 (ICSD #159456)"); //DX20210121 - added two metal-carbides
    vproto.push_back("AB2C_cF32_227_b_c_a;3;8;227;5;1;cF32;a;-;Cs1N2Nb1;Cs1N2Nb1 (ICSD #72546)");

    // -------------------------------------------------------------------------
    // metal prototypes (from DX) //DX20201028
    // -------------------------------------------------------------------------
    // ternaries
    vproto.push_back("ABC2_mC16_8_b_2a_2ab;3;8;8;5;18;mC16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6;-;Bi1Pb1Pd2;Bi1Pb1Pd2 (ICSD #58830)");
    vproto.push_back("A2B4C3_mP9_10_m_2n_am;3;9;10;5;12;mP9;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-;Co2Gd4Mg3;Co2Gd4Mg3 (ICSD #417035)");
    vproto.push_back("A2B5C4_mP11_10_m_a2m_2n;3;11;10;5;14;mP11;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Co2In5Zr4;Co2In5Zr4 (ICSD #55578)");
    vproto.push_back("AB2C2_mP10_11_e_2e_2e;3;10;11;5;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-/-;Ba1Bi2Pd2/Cu1Na2O2;Ba1Bi2Pd2 (ICSD #416299)/Cu1Na2O2 (ICSD #422751)"); //DX20210106 - added metal-oxide
    vproto.push_back("A2BC2_mP10_11_2e_e_2e;3;10;11;5;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Sn2Sr1Zn2;Sn2Sr1Zn2 (ICSD #424108)");
    vproto.push_back("A2B4C_mP14_11_2e_4e_e;3;14;11;5;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Au2In4Yb1;Au2In4Yb1 (ICSD #261042)");
    vproto.push_back("A2BC2_mC10_12_i_a_i;3;5;12;5;8;mC10;a,b/a,c/a,beta,x2,z2,x3,z3;-/-;Cu2Eu1Sn2/Ag2Ni1O2;Cu2Eu1Sn2 (ICSD #182050)/Ag2Ni1O2 (ICSD #172559)"); //DX20210106 - added metal-oxide
    vproto.push_back("ABC_mC12_12_i_i_i;3;6;12;5;10;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-/-;Al1Ce1Co1/K1O1Tl1;Al1Ce1Co1 (ICSD #20632)/K1O1Tl1 (ICSD #1570)"); //DX20210106 - added metal-oxide
    vproto.push_back("A3B2C2_mC14_12_ai_i_i;3;7;12;5;10;mC14;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4;-/-/-/-/-;In3Rh2Sr2/C3Cr2Ho2/C3Ho2Mo2/C3Ce2Mo2/C3Er2Mo2;In3Rh2Sr2 (ICSD #410985)/C3Cr2Ho2 (ICSD #62083)/C3Ho2Mo2 (ICSD #88511)/C3Ce2Mo2 (ICSD #417827)/C3Er2Mo2 (ICSD #617669)"); //DX20210121 - added four metal-carbides
    vproto.push_back("A6B2C_mC18_12_3i_i_a;3;9;12;5;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-;Na6Sn2Zn1;Na6Sn2Zn1 (ICSD #260159)");
    vproto.push_back("A2B5C4_mC22_12_i_a2i_2i;3;11;12;5;14;mC22;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Mn2Sn5Yb4;Mn2Sn5Yb4 (ICSD #419134)");
    vproto.push_back("A3B4C8_mC30_12_ai_hi_2ij;3;15;12;5;16;mC30;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;Ba3Li4Sn8;Ba3Li4Sn8 (ICSD #240016)");
    vproto.push_back("A6B2C7_mC30_12_3i_i_a3i;3;15;12;5;18;mC30;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;Ca6Cu2Sn7;Ca6Cu2Sn7 (ICSD #171243)");
    vproto.push_back("AB2C2_mC20_15_e_f_f;3;10;15;5;11;mC20;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3;-;In1Pd2Sr2;In1Pd2Sr2 (ICSD #391432)");
    vproto.push_back("A2B4C_mC28_15_f_2f_e;3;14;15;5;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-;Na2Sn4Sr1/Li2O4W1/Li2O4W1;Na2Sn4Sr1 (ICSD #261117)/Li2O4W1 (ICSD #1044)/Li2O4W1 (ICSD #14196)"); //DX20210106 - added two metal-oxides
    vproto.push_back("ABC2_oC16_20_b_a_c;3;8;20;5;8;oC16;a,b/a,c/a,x1,y2,x3,y3,z3;-;Ag1Bi1K2;Ag1Bi1K2 (ICSD #1156)");
    vproto.push_back("ABC_oP12_33_a_a_a;3;12;33;5;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ce1Pd1Sn1;Ce1Pd1Sn1 (ICSD #106418)");
    vproto.push_back("ABC_oC12_36_a_a_a;3;6;36;5;9;oC12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-;Ag1Ce1Sn1;Ag1Ce1Sn1 (ICSD #55819)");
    vproto.push_back("ABC2_oC32_36_b_2a_2ab;3;16;36;5;17;oC32;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;-;Bi1Pb1Pd2;Bi1Pb1Pd2 (ICSD #56278)");
    vproto.push_back("ABC_oC18_38_ad_be_bd;3;9;38;5;12;oC18;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6;-;In1Nd1Rh1;In1Nd1Rh1 (ICSD #59430)");
    vproto.push_back("A10BC3_oP14_47_2q2rs_h_bt;3;14;47;5;9;oP14;a,b/a,c/a,z3,z4,z5,z6,z7,z8;-;Ga10Ni1Pr3;Ga10Ni1Pr3 (ICSD #20663)");
    vproto.push_back("ABC2_oP8_51_f_e_be;3;8;51;5;6;oP8;a,b/a,c/a,z2,z3,z4;-;In1La1Ni2;In1La1Ni2 (ICSD #106800)");
    vproto.push_back("A2BC_oP8_51_ae_f_f;3;8;51;5;6;oP8;a,b/a,c/a,z2,z3,z4;-;Co2Ga1La1;Co2Ga1La1 (ICSD #623100)");
    vproto.push_back("A4BC_oP12_51_afj_e_e;3;12;51;5;8;oP12;a,b/a,c/a,z2,z3,z4,x5,z5;-;Al4Co1La1;Al4Co1La1 (ICSD #9986)");
    vproto.push_back("A3B2C_oP12_51_ak_ef_f;3;12;51;5;8;oP12;a,b/a,c/a,z2,z3,z4,y5,z5;-;Au3Rb2Tl1;Au3Rb2Tl1 (ICSD #249924)");
    vproto.push_back("AB4C_oP12_51_e_afj_e;3;12;51;5;8;oP12;a,b/a,c/a,z2,z3,z4,x5,z5;-;Ca1In4Rh1;Ca1In4Rh1 (ICSD #410891)");
    vproto.push_back("ABC_oP12_57_c_d_d;3;12;57;5;8;oP12;a,b/a,c/a,x1,x2,y2,x3,y3;-;Al1Ca1Pd1;Al1Ca1Pd1 (ICSD #370036)");
    vproto.push_back("AB2C_oP16_57_d_2d_c;3;16;57;5;10;oP16;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;-;Bi1K2Sn1;Bi1K2Sn1 (ICSD #107616)");
    vproto.push_back("AB2C_oP8_59_a_e_b;3;8;59;5;7;oP8;a,b/a,c/a,z1,z2,y3,z3;-;Cu1Ni2Ti1;Cu1Ni2Ti1 (ICSD #628580)");
    vproto.push_back("A3B3C_oP14_59_ae_ae_b;3;14;59;5;10;oP14;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5;-;Au3In3Sr1;Au3In3Sr1 (ICSD #245679)");
    vproto.push_back("AB3C4_oP16_59_b_ae_ef;3;16;59;5;11;oP16;a,b/a,c/a,z1,z2,y3,z3,y4,z4,x5,z5;-/-;Au1K3Sn4/Au1Cs3Pb4;Au1K3Sn4 (ICSD #107444)/Au1Cs3Pb4 (ICSD #107448)");
    vproto.push_back("AB2C_oP16_62_c_d_c;3;16;62;5;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;-/-/-;Ca1Cd2Pd1/Ga1Pd2Y1/B1Co2Fe1;Ca1Cd2Pd1 (ICSD #425509)/Ga1Pd2Y1 (ICSD #602790)/B1Co2Fe1 (ICSD #603579)"); //DX20210104 - added ternary-boride
    vproto.push_back("ABC_oC12_63_c_c_c;3;6;63;5;6;oC12;a,b/a,c/a,y1,y2,y3;-/-/-;Mg1Ni1Tb1/Al1B1Mo1/B1Nb1Ni1;Mg1Ni1Tb1 (ICSD #166874)/Al1B1Mo1 (ICSD #16777)/B1Nb1Ni1 (ICSD #43005)"); //DX20210104 - added two ternary borides
    vproto.push_back("ABC_oC12_63_c_c_a;3;6;63;5;5;oC12;a,b/a,c/a,y2,y3;-;Pt1Sc1Zn1;Pt1Sc1Zn1 (ICSD #424600)");
    vproto.push_back("A2BC_oC16_63_f_c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-/E1_{a}/-/-/-;Al2Co1Y1/Al2Cu1Mg1/Ga2Ni1Sc1/Mg2Pd1Sr1/Al2Ce1Pt1;Al2Co1Y1 (ICSD #57645)/Al2Cu1Mg1 (ICSD #57693) and Al2CuMg (part 3)/Ga2Ni1Sc1 (ICSD #103875)/Mg2Pd1Sr1 (ICSD #425493)/Al2Ce1Pt1 (ICSD #658140)"); //DX20210428 - added equivalent part 3 prototype (Al2CuMg, http://aflow.org/prototype-encyclopedia/A2BC_oC16_63_f_c_c.html)
    vproto.push_back("ABC2_oC16_63_a_c_g;3;8;63;5;6;oC16;a,b/a,c/a,y2,x3,y3;-/-/-/-;Au1Bi1Li2/Au1Bi1Li2/Au1Bi1Na2/Au1Bi1K2;Au1Bi1Li2 (ICSD #261785)/Au1Bi1Li2 (ICSD #261786)/Au1Bi1Na2 (ICSD #261788)/Au1Bi1K2 (ICSD #380341)");
    vproto.push_back("ABC2_oC16_63_c_c_2c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,y4;-/-/-;Ba1Cu1Sn2/Co1Ho1Sn2/Ce1Ni1Sn2;Ba1Cu1Sn2 (ICSD #58647)/Co1Ho1Sn2 (ICSD #240098)/Ce1Ni1Sn2 (ICSD #621687)");
    vproto.push_back("AB2C_oC16_63_c_f_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-;La1Mg2Ni1;La1Mg2Ni1 (ICSD #96152)");
    vproto.push_back("A2BC_oC16_63_g_c_c;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,x3,y3;-;Cd2Cu1Er1;Cd2Cu1Er1 (ICSD #99139)");
    vproto.push_back("A2B2C_oC20_63_g_2c_a;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,x4,y4;-;Bi2Cs2Pt1;Bi2Cs2Pt1 (ICSD #658701)");
    vproto.push_back("A4BC_oC24_63_acf_c_c;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,y4,y5,z5;-;Al4Dy1Ni1;Al4Dy1Ni1 (ICSD #57760)");
    vproto.push_back("AB3C2_oC24_63_c_cg_e;3;12;63;5;8;oC24;a,b/a,c/a,y1,y2,x3,x4,y4;-;Ca1Ga3Ni2;Ca1Ga3Ni2 (ICSD #58898)");
    vproto.push_back("A3BC3_oC28_63_cf_a_cf;3;14;63;5;9;oC28;a,b/a,c/a,y2,y3,y4,z4,y5,z5;-/-;Co3Ga1Y3/B3Co1W3;Co3Ga1Y3 (ICSD #10044)/B3Co1W3 (ICSD #25753)"); //DX20210104 - added ternary-boride
    vproto.push_back("ABC2_oC32_63_ac_f_2cf;3;16;63;5;10;oC32;a,b/a,c/a,y2,y3,y4,y5,z5,y6,z6;-;La1Rh1Sn2;La1Rh1Sn2 (ICSD #410732)");
    vproto.push_back("A2B2C_oC10_65_j_i_a;3;5;65;5;5;oC10;a,b/a,c/a,y2,y3;-;Ho2Ni2Pb1;Ho2Ni2Pb1 (ICSD #54612)");
    vproto.push_back("A2B3C_oC12_65_i_cf_a;3;6;65;5;4;oC12;a,b/a,c/a,y4;-;Al2Ni3Pr1;Al2Ni3Pr1 (ICSD #107864)");
    vproto.push_back("A3B2C2_oC14_65_aj_j_i;3;7;65;5;6;oC14;a,b/a,c/a,y2,y3,y4;-;Mg3Ni2Tb2;Mg3Ni2Tb2 (ICSD #240761)");
    vproto.push_back("A2BC_oC16_65_aci_j_i;3;8;65;5;6;oC16;a,b/a,c/a,y3,y4,y5;-;Ga2Nd1Ni1;Ga2Nd1Ni1 (ICSD #103850)");
    vproto.push_back("A2BC6_oC18_65_j_a_2ij;3;9;65;5;7;oC18;a,b/a,c/a,y2,y3,y4,y5;-;Dy2Ni1Sn6;Dy2Ni1Sn6 (ICSD #165783)");
    vproto.push_back("A3B2C7_oC24_65_aj_i_c2ij;3;12;65;5;8;oC24;a,b/a,c/a,y3,y4,y5,y6,y7;-/-;Ce3Ni2Sn7/La3Ni2Sn7;Ce3Ni2Sn7 (ICSD #102239)/La3Ni2Sn7 (ICSD #160913)");
    vproto.push_back("A2B2C_oI10_71_f_h_a;3;5;71;5;5;oI10;a,b/a,c/a,x2,y3;-/-/-;Ce2Ni2Sn1/Ca2Cu2Ga1/K2O2Pd1;Ce2Ni2Sn1 (ICSD #55495)/Ca2Cu2Ga1 (ICSD #58885)/K2O2Pd1 (ICSD #6158)"); //DX20210106 - added metal-oxide
    vproto.push_back("A2B2C3_oI14_71_e_f_af;3;7;71;5;6;oI14;a,b/a,c/a,x2,x3,x4;-;Al2Sn2Sr3;Al2Sn2Sr3 (ICSD #9564)");
    vproto.push_back("A5B2C_oI16_71_an_e_d;3;8;71;5;6;oI16;a,b/a,c/a,x3,x4,y4;-;Al5Ni2Pr1;Al5Ni2Pr1 (ICSD #58046)");
    vproto.push_back("A4B3C4_oI22_71_n_af_eg;3;11;71;5;8;oI22;a,b/a,c/a,x2,x3,y4,x5,y5;-/-;Ag8Ce6Sn8/Ag8Pr6Sn8;Ag8Ce6Sn8 (ICSD #55826)/Ag8Pr6Sn8 (ICSD #55827)");
    vproto.push_back("A4B3C4_oI22_71_n_af_eh;3;11;71;5;8;oI22;a,b/a,c/a,x2,x3,y4,x5,y5;-;Ag4Dy3Sn4;Ag4Dy3Sn4 (ICSD #156968)");
    vproto.push_back("A3B4C4_oI22_71_af_eh_n;3;11;71;5;8;oI22;a,b/a,c/a,x2,x3,y4,x5,y5;-;La3Pd4Zn4;La3Pd4Zn4 (ICSD #182774)");
    vproto.push_back("A6B6C_oI26_71_hk_efg_a;3;13;71;5;7;oI26;a,b/a,c/a,x2,x3,y4,y5;-;Fe6Ga6Sc1;Fe6Ga6Sc1 (ICSD #103465)");
    vproto.push_back("A3B3C8_oI28_71_cg_be_2n;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-/-;Ag3Er3Ga8/Au3Er3Ga8;Ag3Er3Ga8 (ICSD #605107)/Au3Er3Ga8 (ICSD #611824)");
    vproto.push_back("A3B10C_oI28_71_be_g2n_c;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Ba3Hg10In1;Ba3Hg10In1 (ICSD #290305)");
    vproto.push_back("A6B5C3_oI28_71_gn_cn_be;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Cu6Sn5Yb3;Cu6Sn5Yb3 (ICSD #413485)");
    vproto.push_back("A9B3C2_oI28_71_c2n_be_g;3;14;71;5;9;oI28;a,b/a,c/a,x3,y4,x5,y5,x6,y6;-;Ga9Ho3Pd2;Ga9Ho3Pd2 (ICSD #634397)");
    vproto.push_back("A2B2C_oI20_72_j_j_a;3;10;72;5;7;oI20;a,b/a,c/a,x2,y2,x3,y3;-/-;Ba2Bi2Zn1/K2O2Zn1;Ba2Bi2Zn1 (ICSD #421424)/K2O2Zn1 (ICSD #34603)"); //DX20210106 - added metal-oxide
    vproto.push_back("ABC3_tI10_107_a_a_ab;3;5;107;5;6;tI10;a,c/a,z1,z2,z3,z4;-;Ba1Ni1Sn3 (BaNiSn3);Ba1Ni1Sn3 (ICSD #58662) and BaNiSn3 (part 3)"); //DX20210428 - added equivalent part 3 prototype (BaNiSn3, http://aflow.org/prototype-encyclopedia/ABC3_tI10_107_a_a_ab.html)
    vproto.push_back("A3BC_tI10_107_ab_a_a;3;5;107;5;6;tI10;a,c/a,z1,z2,z3,z4;-;Al3Cu1Pr1;Al3Cu1Pr1 (ICSD #290390)");
    vproto.push_back("A3BC_tI10_119_bf_a_c;3;5;119;5;3;tI10;a,c/a,z4;-;In3Mg1Sr1;In3Mg1Sr1 (ICSD #249592)");
    vproto.push_back("AB4C2_tI28_120_d_i_e;3;14;120;5;6;tI28;a,c/a,x2,x3,y3,z3;-;Ce1Ni4Sn2;Ce1Ni4Sn2 (ICSD #102238)");
    vproto.push_back("A4B6C_tI22_121_i_ci_a;3;11;121;5;6;tI22;a,c/a,x3,z3,x4,z4;-;Ru4Sn6Y1;Ru4Sn6Y1 (ICSD #54354)");
    vproto.push_back("AB2C_tP4_123_b_h_a;3;4;123;5;3;tP4;a,c/a,z3;-;Cd1Pt2Zn1;Cd1Pt2Zn1 (ICSD #102057)");
    vproto.push_back("AB4C_tP6_123_b_i_a;3;6;123;5;3;tP6;a,c/a,z3;-;Co1Ga4Hf1;Co1Ga4Hf1 (ICSD #623081)");
    vproto.push_back("AB5C_tP7_123_a_ci_b;3;7;123;5;3;tP7;a,c/a,z4;-;Co1In5Tb1;Co1In5Tb1 (ICSD #623944)");
    vproto.push_back("AB6C_tP8_123_a_hi_b;3;8;123;5;4;tP8;a,c/a,z3,z4;-;Ce1Ga6Pd1;Ce1Ga6Pd1 (ICSD #240161)");
    vproto.push_back("AB8C2_tP11_123_a_ehi_g;3;11;123;5;5;tP11;a,c/a,z3,z4,z5;-/-;Co1Ga8Ho2/Fe1Ga8Ho2;Co1Ga8Ho2 (ICSD #42426)/Fe1Ga8Ho2 (ICSD #180132)");
    vproto.push_back("A2BC2_tP10_127_g_a_h;3;10;127;5;4;tP10;a,c/a,x2,x3;-;Ni2Sn1Zr2;Ni2Sn1Zr2 (ICSD #54303)");
    vproto.push_back("AB2C2_tP10_127_a_g_h;3;10;127;5;4;tP10;a,c/a,x2,x3;-;Al1Cu2Re2;Al1Cu2Re2 (ICSD #57706)");
    vproto.push_back("A2BC_tP8_129_bc_c_a;3;8;129;5;4;tP8;a,c/a,z3,z4;-/-;Sn2Tb1Zn1/Bi2La1Li1;Sn2Tb1Zn1 (ICSD #163425)/Bi2La1Li1 (ICSD #415728)");
    vproto.push_back("AB2C2_tP10_129_c_ac_bc;3;10;129;5;5;tP10;a,c/a,z3,z4,z5;-;Nd1Ni2Sn2;Nd1Ni2Sn2 (ICSD #160053)");
    vproto.push_back("AB2C2_tP10_129_c_bc_ac;3;10;129;5;5;tP10;a,c/a,z3,z4,z5;-;Eu1Sn2Zn2 (Be2CaGe2);Eu1Sn2Zn2 (ICSD #162267) and Be2CaGe2 (part 3)"); //DX20210427 - added equivalent part 3 prototype
    vproto.push_back("A2BC2_tP10_129_bc_c_ac;3;10;129;5;5;tP10;a,c/a,z3,z4,z5;-;Bi2Ce1Ni2;Bi2Ce1Ni2 (ICSD #616562)");
    vproto.push_back("A2BC4_tP14_129_ac_b_j;3;14;129;5;5;tP14;a,c/a,z3,x4,z4;-;Ce2Ru1Zn4;Ce2Ru1Zn4 (ICSD #418547)");
    vproto.push_back("AB2C_tI8_139_b_d_a;3;4;139;5;2;tI8;a,c/a;-;Li1Pd2Tl1;Li1Pd2Tl1 (ICSD #54365)");
    vproto.push_back("A2BC_tI8_139_d_a_b;3;4;139;5;2;tI8;a,c/a;-;Co2Ga1Ni1;Co2Ga1Ni1 (ICSD #169733)");
    vproto.push_back("AB2C2_tI10_139_a_d_e;3;5;139;5;3;tI10;a,c/a,z3;-/-/H5_{9};Ba1Mn2Sn2/C1Li2N2 (CLi2N2)/CaP2U2;Ba1Mn2Sn2 (ICSD #405)/C1Li2N2 (ICSD #200369) and CLi2N2 (part 3)/H59 [Autunite, Ca(UO2)2(PO4)2.10(1/2)H2O] (obsolete)"); //DX20210221 - added carbo-nitride //DX20210427 - added H59 [Autunite, Ca(UO2)2(PO4)2.10(1/2)H2O] (obsolete) //DX20210428 - added equivalent part 3 prototype (CLi2N2, http://aflow.org/prototype-encyclopedia/AB2C2_tI10_139_a_d_e.Li2CN2.html)
    vproto.push_back("A2BC2_tI10_139_d_a_e;3;5;139;5;3;tI10;a,c/a,z3;-/-/-;Ag2Ba1Sn2/Al2Ce1Ga2/Al2Ca1Zn2;Ag2Ba1Sn2 (ICSD #25332)/Al2Ce1Ga2 (ICSD #55789)/Al2Ca1Zn2 (ICSD #57550)");
    vproto.push_back("ABC_tI12_139_c_e_e;3;6;139;5;4;tI12;a,c/a,z2,z3;-;Mg1Pr1Sn1;Mg1Pr1Sn1 (ICSD #182479)");
    vproto.push_back("ABC_tI12_139_e_c_e;3;6;139;5;4;tI12;a,c/a,z2,z3;-;Au1Cu1Sn1;Au1Cu1Sn1 (ICSD #611768)");
    vproto.push_back("A4B2C_tI14_139_h_d_a;3;7;139;5;3;tI14;a,c/a,x3;-;Al4Mo2Yb1;Al4Mo2Yb1 (ICSD #456)");
    vproto.push_back("AB4C2_tI14_139_a_h_d;3;7;139;5;3;tI14;a,c/a,x3;-;Er1Ga4Ti2;Er1Ga4Ti2 (ICSD #630579)");
    vproto.push_back("A5B2C_tI16_139_bg_e_a;3;8;139;5;4;tI16;a,c/a,z3,z4;-;Al5Ni2Zr1;Al5Ni2Zr1 (ICSD #58084)");
    vproto.push_back("A2BC_tI16_139_ce_e_d;3;8;139;5;4;tI16;a,c/a,z3,z4;-;Bi2Sr1Zn1;Bi2Sr1Zn1 (ICSD #41924)");
    vproto.push_back("AB2C_tI16_139_e_ce_d;3;8;139;5;4;tI16;a,c/a,z3,z4;-/-;Ba1Bi2Cd1/Ba1Bi2Zn1;Ba1Bi2Cd1 (ICSD #58635)/Ba1Bi2Zn1 (ICSD #58638)");
    vproto.push_back("A2B5C3_tI20_139_e_bg_ae;3;10;139;5;5;tI20;a,c/a,z3,z4,z5;-;K2Mg5Sn3;K2Mg5Sn3 (ICSD #421342)");
    vproto.push_back("A8BC4_tI26_139_ij_a_f;3;13;139;5;4;tI26;a,c/a,x3,x4;-;Al8Ca1Cu4;Al8Ca1Cu4 (ICSD #57539)");
    vproto.push_back("A4B8C_tI26_139_i_fj_a;3;13;139;5;4;tI26;a,c/a,x3,x4;-;Cr4Fe8Y1;Cr4Fe8Y1 (ICSD #168240)");
    vproto.push_back("A8B4C_tI26_139_ij_f_a;3;13;139;5;4;tI26;a,c/a,x3,x4;-;Al8Mn4Y1;Al8Mn4Y1 (ICSD #57997)");
    vproto.push_back("ABC4_tI24_140_a_c_l;3;12;140;5;4;tI24;a,c/a,x3,z3;-;Ir1Li1Sn4;Ir1Li1Sn4 (ICSD #172149)");
    vproto.push_back("A2BC4_tI28_140_h_a_k;3;14;140;5;5;tI28;a,c/a,x2,x3,y3;-;Bi2Mn1Ti4 (Sb2SiV4);Bi2Mn1Ti4 (ICSD #150145) and Sb2SiV4 (part 3)"); //DX20210427 - added equivalent part 3 prototype
    vproto.push_back("A4B2C_tI28_140_l_h_a;3;14;140;5;5;tI28;a,c/a,x2,x3,z3;-;Cu4Sn2Sr1;Cu4Sn2Sr1 (ICSD #182104)");
    vproto.push_back("AB5C2_tI32_140_a_cl_h;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Bi1Er5Pt2;Bi1Er5Pt2 (ICSD #107217)");
    vproto.push_back("A3B4C_tI32_140_ah_l_c;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Bi3In4Pb1;Bi3In4Pb1 (ICSD #616725)");
    vproto.push_back("AB5C2_tI32_140_a_bk_h;3;16;140;5;5;tI32;a,c/a,x3,x4,y4;-;Ga1Nb5Sn2;Ga1Nb5Sn2 (ICSD #103842)");
    vproto.push_back("A2BC5_tI32_140_h_a_cl;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Au2Bi1Tb5;Au2Bi1Tb5 (ICSD #156956)");
    vproto.push_back("AB2C5_tI32_140_a_h_cl;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;-;Bi1Co2Ho5;Bi1Co2Ho5 (ICSD #161658)");
    vproto.push_back("A4B7C3_hR14_146_ab_a2b_b;3;14;146;5;16;hR14;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Ga4Pd7Zn3;Ga4Pd7Zn3 (ICSD #103911)");
    vproto.push_back("ABC_hP9_156_3a_b2c_2bc;3;9;156;5;11;hP9;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8,z9;-;Ca1Li1Sn1;Ca1Li1Sn1 (ICSD #58911)");
    vproto.push_back("A5BC4_hR10_160_5a_a_4a;3;10;160;5;12;hR10;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10;-;Li5Na1Sn4;Li5Na1Sn4 (ICSD #12142)");
    vproto.push_back("A2BC_hR4_166_c_b_a;3;4;166;5;3;hR4;a,c/a,x3;-;Cu2Ga1Sr1;Cu2Ga1Sr1 (ICSD #102938)");
    vproto.push_back("AB2C4_hR7_166_a_c_2c;3;7;166;5;5;hR7;a,c/a,x2,x3,x4;-;Au1Ni2Sn4;Au1Ni2Sn4 (ICSD #150127)");
    vproto.push_back("AB4C6_hR11_166_a_2c_h;3;11;166;5;6;hR11;a,c/a,x2,x3,x4,z4;-;Gd1Sn4Ti6;Gd1Sn4Ti6 (ICSD #166193)");
    vproto.push_back("A2BC_hR12_166_h_bc_ac;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Al2Cu1Yb1;Al2Cu1Yb1 (ICSD #604213)");
    vproto.push_back("A7B3C2_hR12_166_bh_ac_c;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Al7Ca3Cu2;Al7Ca3Cu2 (ICSD #57538)");
    vproto.push_back("AB2C9_hR12_166_a_c_bch;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;La1Mg2Ni9;La1Mg2Ni9 (ICSD #55614)");
    vproto.push_back("A9B2C_hR12_166_eh_c_a;3;12;166;5;5;hR12;a,c/a,x2,x4,z4;-;Ni3Ti0.67Zr0.33;Ni3Ti0.67Zr0.33 (ICSD #646978)");
    vproto.push_back("AB2C_hR12_166_bc_h_ac;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Ag1Al2Pr1;Ag1Al2Pr1 (ICSD #604688)");
    vproto.push_back("A2B7C3_hR12_166_c_bh_ac;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;Ag2Al7Ca3;Ag2Al7Ca3 (ICSD #104173)");
    vproto.push_back("A7B4C2_hR13_166_ah_2c_c;3;13;166;5;7;hR13;a,c/a,x2,x3,x4,x5,z5;-;Au7Rb4Sn2;Au7Rb4Sn2 (ICSD #58581)");
    vproto.push_back("A3B2C2_hR14_166_3c_2c_abc;3;14;166;5;8;hR14;a,c/a,x3,x4,x5,x6,x7,x8;-;Au3Sn2Yb2;Au3Sn2Yb2 (ICSD #710044)");
    vproto.push_back("ABC_hP6_186_b_a_b;3;6;186;5;5;hP6;a,c/a,z1,z2,z3;-/-;Au1Pr1Sn1/Bi1Li1Zn1;Au1Pr1Sn1 (ICSD #54997)/Bi1Li1Zn1 (ICSD #100115)");
    vproto.push_back("ABC_hP6_186_a_b_b;3;6;186;5;5;hP6;a,c/a,z1,z2,z3;-;Ce1Cu1Sn1;Ce1Cu1Sn1 (ICSD #156394)");
    vproto.push_back("A4BC_hP12_186_ac_b_b;3;12;186;5;7;hP12;a,c/a,z1,z2,z3,x4,z4;-;Cu4In1Mn1;Cu4In1Mn1 (ICSD #424277)");
    vproto.push_back("A3BC4_hP16_186_a2b_b_2a2b;3;16;186;5;10;hP16;a,c/a,z1,z2,z3,z4,z5,z6,z7,z8;-;Au3Li1Sn4;Au3Li1Sn4 (ICSD #412207)");
    vproto.push_back("ABC_hP6_187_i_h_ab;3;6;187;5;4;hP6;a,c/a,z3,z4;-;Au1Ga1Zr1;Au1Ga1Zr1 (ICSD #156264)");
    vproto.push_back("A4BC_hP6_187_hi_b_a;3;6;187;5;4;hP6;a,c/a,z3,z4;-;Ga4Li1Y1;Ga4Li1Y1 (ICSD #98666)");
    vproto.push_back("ABC_hP9_187_dh_ai_eg;3;9;187;5;5;hP9;a,c/a,z4,z5,z6;-;Ca1Li1Pb1;Ca1Li1Pb1 (ICSD #409533)");
    vproto.push_back("A5B5C_hP11_187_hk_cgi_a;3;11;187;5;6;hP11;a,c/a,z3,z4,z5,x6;-;Al5Ba5Pb1;Al5Ba5Pb1 (ICSD #416337)");
    vproto.push_back("A4B6C_hP11_187_ck_jk_a;3;11;187;5;5;hP11;a,c/a,x3,x4,x5;-;Au4In6K1;Au4In6K1 (ICSD #249520)");
    vproto.push_back("A3B2C3_hP8_189_f_c_g;3;8;189;5;4;hP8;a,c/a,x2,x3;-;In3Rh2Ti3;In3Rh2Ti3 (ICSD #410967)");
    vproto.push_back("ABC_hP9_189_g_f_ad;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Bi1Dy1Rh1 (AlNiZr);Bi1Dy1Rh1 (ICSD #51845) and AlNiZr (part 3)"); //DX20210428 - added equivalent part 3 prototype (AlNiZr, http://aflow.org/prototype-encyclopedia/ABC_hP9_189_g_ad_f.html)
    vproto.push_back("AB2C6_hP9_189_b_c_fg;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Co1Ga2Zr6;Co1Ga2Zr6 (ICSD #20876)");
    vproto.push_back("A5B3C_hP9_189_cf_g_b;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Ga5Lu3Ni1;Ga5Lu3Ni1 (ICSD #634575)");
    vproto.push_back("ABC_hP9_189_g_bc_f;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Pb1Pd1Y1;Pb1Pd1Y1 (ICSD #54315)");
    vproto.push_back("A2BC6_hP9_189_c_b_fg;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Bi2Fe1Ho6;Bi2Fe1Ho6 (ICSD #96253)");
    vproto.push_back("ABC_hP9_189_f_bc_g;3;9;189;5;4;hP9;a,c/a,x3,x4;-/-;Mg1Sn1Yb1/Ga1Ni1Ti1;Mg1Sn1Yb1 (ICSD #54344)/Ga1Ni1Ti1 (ICSD #103885)");
    vproto.push_back("ABC_hP9_189_f_g_bc;3;9;189;5;4;hP9;a,c/a,x3,x4;-/-/-;Al1Dy1Ni1/Al1Hf1Pt1/Al1Nd1Ni1;Al1Dy1Ni1 (ICSD #107416)/Al1Hf1Pt1 (ICSD #608141)/Al1Nd1Ni1 (ICSD #608770)");
    vproto.push_back("ABC_hP9_189_bc_f_g;3;9;189;5;4;hP9;a,c/a,x3,x4;-;Cu1In1Tb1;Cu1In1Tb1 (ICSD #628137)");
    vproto.push_back("A3B2C_hP6_191_g_c_a;3;6;191;5;2;hP6;a,c/a;-/-;Ag3Al2La1/Co3Ga2Tb1;Ag3Al2La1 (ICSD #57329)/Co3Ga2Tb1 (ICSD #102452)");
    vproto.push_back("A2BC3_hP6_191_c_a_g;3;6;191;5;2;hP6;a,c/a;-/-;Cu2Dy1In3/Al2Ce1Pt3;Cu2Dy1In3 (ICSD #108369)/Al2Ce1Pt3 (ICSD #658142)");
    vproto.push_back("A9BC2_hP12_191_fm_a_c;3;12;191;5;3;hP12;a,c/a,x4;-;Al9Ba1Fe2;Al9Ba1Fe2 (ICSD #57518)");
    vproto.push_back("A6BC6_hP13_191_i_a_cde;3;13;191;5;4;hP13;a,c/a,z4,z5;-;Mn6Sc1Sn6;Mn6Sc1Sn6 (ICSD #54273)");
    vproto.push_back("A6B6C_hP13_191_i_cde_a;3;13;191;5;4;hP13;a,c/a,z4,z5;-;Mn6Sn6Tb1;Mn6Sn6Tb1 (ICSD #54277)");
    vproto.push_back("A3B3C2_hP16_193_g_g_d;3;16;193;5;4;hP16;a,c/a,x2,x3;-;Ga3Hf3Nb2;Ga3Hf3Nb2 (ICSD #103733)");
    vproto.push_back("ABC_hP6_194_d_a_c;3;6;194;5;2;hP6;a,c/a;-;Pb1Sr1Zn1;Pb1Sr1Zn1 (ICSD #54319)");
    vproto.push_back("ABC_hP6_194_c_a_d;3;6;194;5;2;hP6;a,c/a;-/-;Bi1Ca1Cu1/Cu1Gd1Sn1;Bi1Ca1Cu1 (ICSD #57018)/Cu1Gd1Sn1 (ICSD #600994)");
    vproto.push_back("ABC_hP6_194_a_c_d;3;6;194;5;2;hP6;a,c/a;-;Ba1Pb1Zn1;Ba1Pb1Zn1 (ICSD #106315)");
    vproto.push_back("ABC_hP6_194_a_d_c;3;6;194;5;2;hP6;a,c/a;-;Ca1Hg1Pb1;Ca1Hg1Pb1 (ICSD #106345)");
    vproto.push_back("AB2C_hP8_194_b_f_c;3;8;194;5;3;hP8;a,c/a,z3;-;Al1Pt2Zr1;Al1Pt2Zr1 (ICSD #58139)");
    vproto.push_back("AB2C_hP8_194_a_f_c;3;8;194;5;3;hP8;a,c/a,z3;-;In1Pt2Y1;In1Pt2Y1 (ICSD #59505)");
    vproto.push_back("A2BC_hP8_194_f_a_c;3;8;194;5;3;hP8;a,c/a,z3;-;Cu2Li1Sn1;Cu2Li1Sn1 (ICSD #150602)");
    vproto.push_back("A3BC2_hP12_194_h_a_f;3;12;194;5;4;hP12;a,c/a,z2,x3;-;Al3Ru1Sc2;Al3Ru1Sc2 (ICSD #58159)");
    vproto.push_back("ABC_hP12_194_f_f_ab;3;12;194;5;4;hP12;a,c/a,z3,z4;-;Ga1Sn1Sr1;Ga1Sn1Sr1 (ICSD #66003)");
    vproto.push_back("ABC_hP12_194_ab_f_f;3;12;194;5;4;hP12;a,c/a,z3,z4;-/-;Gd1Sn1Zn1/La1Sn1Zn1;Gd1Sn1Zn1 (ICSD #152621)/La1Sn1Zn1 (ICSD #152623)");
    vproto.push_back("A3B3C2_hP16_194_h_af_f;3;16;194;5;5;hP16;a,c/a,z2,z3,x4;-;Al2Ba2Ga1.33;Al2Ba2Ga1.33 (ICSD #415931)");
    vproto.push_back("AB2C_cP16_215_e_bce_ad;3;16;215;5;3;cP16;a,x5,x6;-;Ir1Li2Mg1;Ir1Li2Mg1 (ICSD #180119)");
    vproto.push_back("ABC2_cF16_216_a_d_bc;3;4;216;5;1;cF16;a;-;Ag1Al1Li2 (CuHg2Ti);Ag1Al1Li2 (ICSD #57330) and CuHg2Ti Inverse Heusler (part )"); //DX20210428 - added equivalent part 3 prototype (CuHg2Ti, inverse Heusler, http://aflow.org/prototype-encyclopedia/AB2C_cF16_216_b_ad_c.html)
    vproto.push_back("AB4C_cF24_216_a_e_c;3;6;216;5;2;cF24;a,x3;-;In1Ni4Zr1;In1Ni4Zr1 (ICSD #59462)");
    vproto.push_back("A3B4C4_cI22_217_b_c_c;3;11;217;5;3;cI22;a,x2,x3;-;Ca3Ga4Ni4;Ca3Ga4Ni4 (ICSD #58899)");
    vproto.push_back("A3BC_cP15_221_ag_c_d;3;15;221;5;2;cP15;a,x4;-;Ga3Os1Tb1;Ga3Os1Tb1 (ICSD #412144)");

    // -------------------------------------------------------------------------
    // metal-boride prototypes (from DX) //DX20210104
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A4B_mC10_12_j_a;2;5;12;4;7;mC10;a,b/a,c/a,beta,x2,y2,z2;-;B4Mn1;B4Mn1 (ICSD #15079)");
    vproto.push_back("A3B4_mC28_15_ef_2f;2;14;15;4;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-;B3Ni4/B3Ni4;B3Ni4 (ICSD #24308)/B3Ni4 (ICSD #150561)");
    vproto.push_back("A4B_oP10_58_2g_a;2;10;58;4;7;oP10;a,b/a,c/a,x2,y2,x3,y3;-;B4Cr1;B4Cr1 (ICSD #186851)");
    vproto.push_back("A2B_oP6_59_e_a;2;6;59;4;6;oP6;a,b/a,c/a,z1,y2,z2;-;B2Ru1 (B2Ru);B2Ru1 (ICSD #31871) and B2Ru (part 3)"); //DX20210428 - added equivalent part 3 prototype (B2Ru, http://aflow.org/prototype-encyclopedia/A2B_oP6_59_f_a.html)
    vproto.push_back("A2B_oP12_62_d_c;2;12;62;4;8;oP12;a,b/a,c/a,x1,z1,x2,y2,z2;-;B2Fe1;B2Fe1 (ICSD #425310)");
    vproto.push_back("AB_oC8_63_a_c;2;4;63;4;4;oC8;a,b/a,c/a,y2;-;B1Rh1;B1Rh1 (ICSD #150732)");
    vproto.push_back("AB3_oC16_63_c_cf;2;8;63;4;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-;B1Re3 (BRe3);B1Re3 (ICSD #43662) and Re3B (part 3)"); //DX20210428 - added equivalent part 3 prototype (Re3B, http://aflow.org/prototype-encyclopedia/AB3_oC16_63_c_cf.html)
    vproto.push_back("A3B2_oC20_63_3c_2c;2;10;63;4;8;oC20;a,b/a,c/a,y1,y2,y3,y4,y5;-;B3V2;B3V2 (ICSD #79258)");
    vproto.push_back("A6B5_oC22_65_i2j_aij;2;11;65;4;8;oC22;a,b/a,c/a,y2,y3,y4,y5,y6;-;B6Ta5;B6Ta5 (ICSD #68538)");
    vproto.push_back("A4B_oI10_71_n_a;2;5;71;4;5;oI10;a,b/a,c/a,x2,y2;-;B4Cr1;B4Cr1 (ICSD #24353)");
    vproto.push_back("A4B3_oI14_71_ef_af;2;7;71;4;6;oI14;a,b/a,c/a,x2,x3,x4;-/-;B4Nb3/B4Mn3;B4Nb3 (ICSD #76631)/B4Mn3 (ICSD #614732)");
    vproto.push_back("A11B_tI24_119_aegi_d;2;12;119;4;6;tI24;a,c/a,z3,x4,x5,z5;-/-;B11Li1/B11Li1;B11Li1 (ICSD #164844)/B11Li1 (ICSD #164845)");
    vproto.push_back("AB2_tI12_121_ab_i;2;6;121;4;4;tI12;a,c/a,x3,z3;C17;B1Fe2 (BFe2);B1Fe2 (ICSD #16809) Fe2B (C17, obsolete, part 3)"); //DX20210428 - added equivalent part 3 prototype (Fe2B, C17, obsolete, http://aflow.org/prototype-encyclopedia/AB2_tI12_121_ab_i.html)
    vproto.push_back("A12B_tI26_139_jm_a;2;13;139;4;5;tI26;a,c/a,x2,x3,z3;-;B12Sc1;B12Sc1 (ICSD #615424)");
    vproto.push_back("AB2_tI12_140_b_h;2;6;140;4;3;tI12;a,c/a,x2;-/-;B1Fe2/B1Fe2;B1Fe2 (ICSD #160789)/B1Fe2 (ICSD #160790)");
    vproto.push_back("AB2_tI12_140_c_h;2;6;140;4;3;tI12;a,c/a,x2;-;B1Ti2;B1Ti2 (ICSD #189385)");
    vproto.push_back("A4B5_hR9_160_ab_2ab;2;9;160;4;9;hR9;a,c/a,x1,x2,x3,x4,z4,x5,z5;-;B4Li5;B4Li5 (ICSD #23041)");
    vproto.push_back("A11B_hR12_160_2a3b_a;2;12;160;4;11;hR12;a,c/a,x1,x2,x3,x4,z4,x5,z5,x6,z6;-;B11Li1;B11Li1 (ICSD #164842)");
    vproto.push_back("A2B_hR6_166_2c_c;2;6;166;4;5;hR6;a,c/a,x1,x2,x3;-;B4Mo2;B4Mo2 (ICSD #39554)");
    vproto.push_back("A3B_hR8_166_f_c;2;8;166;4;4;hR8;a,c/a,x1,x2;-;B3Mo1;B3Mo1 (ICSD #167734)");
    vproto.push_back("A2B_hP3_191_d_a;2;3;191;4;2;hP3;a,c/a;-/-;B2Zr1/B2Ca1/B2Nb1;B2Zr1 (ICSD #169458)/B2Ca1 (ICSD #186764)/B2Nb1 (ICSD #614887)");
    vproto.push_back("A2B_hP6_194_f_c;2;6;194;4;3;hP6;a,c/a,z2;-/-;B2Re1/C2Re1;B2Re1 (ICSD #23871)/C2Re1 (ICSD #184660)"); //DX20210120 - added metal-carbide
    vproto.push_back("AB_hP8_194_bc_f;2;8;194;4;3;hP8;a,c/a,z3;-;B2Li2;B2Li2 (ICSD #2)");
    vproto.push_back("A3B_hP8_194_af_c;2;8;194;4;3;hP8;a,c/a,z3;-;B3Re1 (B3Re);B3Re1 (ICSD #24361) and B3Re (part 3)"); //DX20210428 - added equivalent part 3 prototype (B3Re, http://aflow.org/prototype-encyclopedia/A3B_hP8_194_af_c.html)
    vproto.push_back("A4B_hP10_194_2f_c;2;10;194;4;4;hP10;a,c/a,z2,z3;-;B4Mo1;B4Mo1 (ICSD #167735)");
    vproto.push_back("A3B2_hP10_194_bf_f;2;10;194;4;4;hP10;a,c/a,z2,z3;-;B3Ru2;B3Ru2 (ICSD #108082)");
    vproto.push_back("A2B_hP12_194_bcf_f;2;12;194;4;4;hP12;a,c/a,z3,z4;-;B2W1;B2W1 (ICSD #23716)");
    vproto.push_back("A5B2_hP14_194_fg_e;2;14;194;4;4;hP14;a,c/a,z1,z2;-;B5Mo2;B5Mo2 (ICSD #167733)");
    vproto.push_back("A11B_cF48_216_aeg_b;2;12;216;4;3;cF48;a,x3,x4;-;B11Li1;B11Li1 (ICSD #164841)");
    vproto.push_back("A11B_cF48_216_aef_d;2;12;216;4;3;cF48;a,x3,x4;-;B11Li1;B11Li1 (ICSD #164843)");
    // ternaries
    vproto.push_back("A4B2C_aP7_2_2i_i_a;3;7;2;5;15;aP7;a,b/a,c/a,alpha,beta,gamma,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;B4Co2Mo1;B4Co2Mo1 (ICSD #44175)");
    vproto.push_back("A5B2C3_mC20_5_b2c_c_ac;3;10;5;5;18;mC20;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;B5Ca2Os3;B5Ca2Os3 (ICSD #59227)");
    vproto.push_back("A3BC_mP10_11_3e_e_e;3;10;11;5;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;B3Er1Mo1;B3Er1Mo1 (ICSD #65932)");
    vproto.push_back("A2BC3_mC12_12_g_a_df;3;6;12;5;5;mC12;a,b/a,c/a,beta,y4;-;B2Er1Ir3;B2Er1Ir3 (ICSD #44230)");
    vproto.push_back("A4B2C5_mC22_12_2i_i_aj;3;11;12;5;13;mC22;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,y5,z5;-/-;B4La2Ni5/B4La2Ni5;B4La2Ni5 (ICSD #63501)/B4La2Ni5 (ICSD #170618)");
    vproto.push_back("A2B2C_mC20_15_f_f_e;3;10;15;5;11;mC20;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3;-;B2Ni2Tb1;B2Ni2Tb1 (ICSD #57003)");
    vproto.push_back("ABC_oF48_43_b_b_b;3;12;43;5;12;oF48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;B1Cu1Ir1;B1Cu1Ir1 (ICSD #75029)");
    vproto.push_back("AB4C_oP6_47_d_rt_a;3;6;47;5;5;oP6;a,b/a,c/a,z3,z4;-;Al1B4Yb1;Al1B4Yb1 (ICSD #181368)");
    vproto.push_back("AB4C3_oP8_47_a_qr_dt;3;8;47;5;6;oP8;a,b/a,c/a,z3,z4,z5;-;Al1B4Cr3;Al1B4Cr3 (ICSD #20082)");
    vproto.push_back("A3B4C_oP16_47_qrs_qrst_ad;3;16;47;5;10;oP16;a,b/a,c/a,z3,z4,z5,z6,z7,z8,z9;-;B3Ir4Zn1;B3Ir4Zn1 (ICSD #71640)");
    vproto.push_back("A2B2C3_oP14_55_h_h_ag;3;14;55;5;9;oP14;a,b/a,c/a,x2,y2,x3,y3,x4,y4;-;B2Li2Rh3;B2Li2Rh3 (ICSD #417442)");
    vproto.push_back("A2BC_oP16_62_2c_c_c;3;16;62;5;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;-/-;B2Co1Nb1/C2Dy1Mo1;B2Co1Nb1 (ICSD #41893)/C2Dy1Mo1 (ICSD #617606)"); //DX20210121 - added metal-carbide
    vproto.push_back("A2BC_oP16_62_d_c_c;3;16;62;5;10;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,y3,z3;-/-/-;B2Ce1Ru1/C2Cs1K1/C2Cs1Rb1;B2Ce1Ru1 (ICSD #612853)/C2Cs1K1 (ICSD #189822)/C2Cs1Rb1 (ICSD #189823)"); //DX20210121 - added two metal-carbides
    vproto.push_back("A3BC_oC20_63_3c_c_c;3;10;63;5;8;oC20;a,b/a,c/a,y1,y2,y3,y4,y5;-;B3Co1V1;B3Co1V1 (ICSD #44188)");
    vproto.push_back("A4BC3_oC32_63_2cf_a_cf;3;16;63;5;10;oC32;a,b/a,c/a,y2,y3,y4,y5,z5,y6,z6;-;B4Mg1Os3;B4Mg1Os3 (ICSD #91243)");
    vproto.push_back("AB2C2_oC10_65_a_i_j;3;5;65;5;5;oC10;a,b/a,c/a,y2,y3;-;Al1B2Cr2;Al1B2Cr2 (ICSD #20083)");
    vproto.push_back("A2BC3_oC12_65_i_a_cf;3;6;65;5;4;oC12;a,b/a,c/a,y4;-;B2Ce1Rh3;B2Ce1Rh3 (ICSD #40777)");
    vproto.push_back("A2B2C3_oC14_65_j_j_ai;3;7;65;5;6;oC14;a,b/a,c/a,y2,y3,y4;-;Al2B2Ru3;Al2B2Ru3 (ICSD #43843)");
    vproto.push_back("A6B2C3_oC22_65_gp_j_ah;3;11;65;5;8;oC22;a,b/a,c/a,x2,x3,y4,x5,y5;-;B6Er2Ni3;B6Er2Ni3 (ICSD #65737)");
    vproto.push_back("ABC4_oC24_65_i_bc_jkm;3;12;65;5;7;oC24;a,b/a,c/a,y3,y4,z5,z6;-;B1Nd1Ni4;B1Nd1Ni4 (ICSD #20882)");
    vproto.push_back("AB4C_oC24_65_h_gip_j;3;12;65;5;9;oC24;a,b/a,c/a,x1,x2,y3,y4,x5,y5;-;Al1B4Yb1;Al1B4Yb1 (ICSD #167533)");
    vproto.push_back("A4B2C5_oF44_69_o_g_acf;3;11;69;5;6;oF44;a,b/a,c/a,x4,x5,y5;-;B4Ca2Rh5;B4Ca2Rh5 (ICSD #66768)");
    vproto.push_back("A2B2C_oF40_70_g_e_a;3;10;70;5;5;oF40;a,b/a,c/a,x2,z3;-;B2Rh2Sr1;B2Rh2Sr1 (ICSD #8153)");
    vproto.push_back("A2BC2_oF40_70_f_a_e;3;10;70;5;5;oF40;a,b/a,c/a,x2,y3;-;B2Nd1Os2;B2Nd1Os2 (ICSD #601369)");
    vproto.push_back("ABC_oF48_70_e_f_c;3;12;70;5;5;oF48;a,b/a,c/a,x2,y3;-;B1Ir1Li1;B1Ir1Li1 (ICSD #75027)");
    vproto.push_back("A2BC2_oI10_71_h_a_f;3;5;71;5;5;oI10;a,b/a,c/a,x2,y3;-;B2Co1W2;B2Co1W2 (ICSD #16776)");
    vproto.push_back("A4B2C_oI14_71_ef_f_a;3;7;71;5;6;oI14;a,b/a,c/a,x2,x3,x4;-/-;B4Fe2Mo1/B4Mn2W1;B4Fe2Mo1 (ICSD #44293)/B4Mn2W1 (ICSD #614782)");
    vproto.push_back("A6BC2_oI18_71_hn_a_f;3;9;71;5;7;oI18;a,b/a,c/a,x2,y3,x4,y4;-/-;B6Ce1Cr2/B6Ce1Cr2;B6Ce1Cr2 (ICSD #16203)/B6Ce1Cr2 (ICSD #81544)");
    vproto.push_back("ABC3_tP5_99_a_b_ac;3;5;99;5;6;tP5;a,c/a,z1,z2,z3,z4;-;B1Ce1Pt3;B1Ce1Pt3 (ICSD #95049)");
    vproto.push_back("A2BC2_tI10_139_e_a_d;3;5;139;5;3;tI10;a,c/a,z3;-;B2Ba1Rh2;B2Ba1Rh2 (ICSD #8155)");
    vproto.push_back("A2B2C_tI10_139_d_e_a;3;5;139;5;3;tI10;a,c/a,z3;-;B2Ir2Zn1;B2Ir2Zn1 (ICSD #71642)");
    vproto.push_back("A2B2C_tI10_139_e_d_a;3;5;139;5;3;tI10;a,c/a,z3;-;B2Co2Dy1 (Cr2Si2Th);B2Co2Dy1 (ICSD #612907) and Cr2Si2Th (part 3)"); //DX20210428 - added equivalent part 3 prototype (Cr2Si2Th, http://aflow.org/prototype-encyclopedia/A2B2C_tI10_139_d_e_a.html)
    vproto.push_back("A2B3C7_hR12_166_c_bc_ah;3;12;166;5;6;hR12;a,c/a,x3,x4,x5,z5;-;B2Ca3Ni7;B2Ca3Ni7 (ICSD #36505)");
    vproto.push_back("A6B3C4_hR13_166_g_ac_2c;3;13;166;5;6;hR13;a,c/a,x2,x3,x4,x5;-;B6Fe3Pr4;B6Fe3Pr4 (ICSD #614152)");
    vproto.push_back("A4BC3_hP16_176_bh_c_h;3;16;176;5;6;hP16;a,c/a,x3,y3,x4,y4;-;B4Hf1Ir3;B4Hf1Ir3 (ICSD #1518)");
    vproto.push_back("ABC2_hP12_180_d_c_i;3;12;180;5;3;hP12;a,c/a,x3;-;B1Ce1Pt2;B1Ce1Pt2 (ICSD #90403)");
    vproto.push_back("ABC_hP9_189_ac_f_g;3;9;189;5;4;hP9;a,c/a,x3,x4;-;B1Fe1Nb1;B1Fe1Nb1 (ICSD #20298)");
    vproto.push_back("A2B6C5_hP13_189_ab_i_df;3;13;189;5;5;hP13;a,c/a,x4,x5,z5;-;B2Rh6Sn5;B2Rh6Sn5 (ICSD #77352)");
    vproto.push_back("ABC3_hP15_189_bc_f_fk;3;15;189;5;6;hP15;a,c/a,x3,x4,x5,y5;-;B1Li1Pt3;B1Li1Pt3 (ICSD #68091)");
    vproto.push_back("ABC3_hP10_191_c_ab_i;3;10;191;5;3;hP10;a,c/a,z4;-;B1Na1Pt3;B1Na1Pt3 (ICSD #68092)");
    vproto.push_back("AB4C_hP12_191_d_ci_ab;3;12;191;5;3;hP12;a,c/a,z5;-;B1Ni4Y1;B1Ni4Y1 (ICSD #16515)");
    vproto.push_back("ABC4_hP12_191_c_ab_di;3;12;191;5;3;hP12;a,c/a,z5;-;B1Ca1Ni4;B1Ca1Ni4 (ICSD #36504)");
    vproto.push_back("A2B7C_cF40_225_c_bd_a;3;10;225;5;1;cF40;a;-;B2Pd7Y1;B2Pd7Y1 (ICSD #260147)");
    // -------------------------------------------------------------------------
    // metal-oxide prototypes (from DX) //DX 20210106
    // -------------------------------------------------------------------------
    // binaries
    // -------------------------------------------------------------------------
    vproto.push_back("AB2_mP12_4_2a_4a;2;12;4;4;22;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Mo1O2;Mo1O2 (ICSD #36263)");
    vproto.push_back("A4B_mC10_5_2c_a;2;5;5;4;11;mC10;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3;-;O4Os1;O4Os1 (ICSD #24672)");
    vproto.push_back("A2B5_mC14_5_c_a2c;2;7;5;4;14;mC14;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Nb2O5;Nb2O5 (ICSD #51176)");
    vproto.push_back("AB2_mC6_8_a_2a;2;3;8;4;10;mC6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;Co1O2;Co1O2 (ICSD #95440)");
    vproto.push_back("AB_mC8_9_a_a;2;4;9;4;10;mC8;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2;-;Cu1O1;Cu1O1 (ICSD #69757)");
    vproto.push_back("A2B_mP6_10_mn_ce;2;6;10;4;8;mP6;a,b/a,c/a,beta,x3,z3,x4,z4;-;O2V1;O2V1 (ICSD #89471)");
    vproto.push_back("A2B_mP12_10_2mo_im;2;12;10;4;14;mP12;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;O2V1;O2V1 (ICSD #1503)");
    vproto.push_back("AB3_mP8_11_e_3e;2;8;11;4;12;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;Mo1O3;Mo1O3 (ICSD #80577)");
    vproto.push_back("A3B2_mP10_11_3e_2e;2;10;11;4;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;O3Pb2;O3Pb2 (ICSD #36243)");
    vproto.push_back("A2B_mP12_11_4e_2e;2;12;11;4;16;mP12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;O2Ti1;O2Ti1 (ICSD #657748)");
    vproto.push_back("A3B4_mP14_11_3e_4e;2;14;11;4;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;O3Tl4;O3Tl4 (ICSD #23478)");
    vproto.push_back("A5B2_mP14_11_5e_2e;2;14;11;4;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;O5V2;O5V2 (ICSD #59960)");
    vproto.push_back("AB2_mC12_12_i_2i;2;6;12;4;10;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-/-;Nd2O4/O1Tl2;Nd2O4 (ICSD #62228)/O1Tl2 (ICSD #77699)");
    vproto.push_back("AB2_mC18_12_ai_3i;2;9;12;4;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-/-;Mn1O2/NbTe2;Mn1O2 (ICSD #150462)/NbTe2"); //DX20210427 - added NbTe2 from part 3
    vproto.push_back("A3B7_mC20_12_di_a3i;2;10;12;4;12;mC20;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6;-;O3V7;O3V7 (ICSD #77706)");
    vproto.push_back("A3B7_mC20_12_ai_d3i;2;10;12;4;12;mC20;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6;-;O6V14;O6V14 (ICSD #163494)");
    vproto.push_back("A3B8_mC22_12_ai_2ij;2;11;12;4;13;mC22;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;Cr3O8;Cr3O8 (ICSD #155847)");
    vproto.push_back("A2B_mC24_12_4i_2i;2;12;12;4;16;mC24;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-/-;O2V1/O2Ti1/O2V1;O2V1 (ICSD #199)/O2Ti1 (ICSD #57154)/O2V1 (ICSD #73855)");
    vproto.push_back("A2B_mC24_12_2ij_gi;2;12;12;4;14;mC24;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;O2V1;O2V1 (ICSD #34416)");
    vproto.push_back("A5B8_mC26_12_ahi_2ij;2;13;12;4;14;mC26;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;Mn5O8;Mn5O8 (ICSD #16956)");
    vproto.push_back("A4B3_mP14_13_2g_ag;2;14;13;4;13;mP14;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;O4Sn3;O4Sn3 (ICSD #174299)");
    vproto.push_back("AB_mP8_14_ad_e;2;8;14;4;7;mP8;a,b/a,c/a,beta,x3,y3,z3;-;Ag1O1;Ag1O1 (ICSD #43741)");
    vproto.push_back("A3B4_mP14_14_ae_2e;2;14;14;4;13;mP14;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Ag3O4;Ag3O4 (ICSD #59225)");
    vproto.push_back("AB2_mC12_15_c_f;2;6;15;4;7;mC12;a,b/a,c/a,beta,x2,y2,z2;-;K1O2;K1O2 (ICSD #38230)");
    vproto.push_back("AB2_mC24_15_ae_2f;2;12;15;4;11;mC24;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;-;Bi2O4;Bi2O4 (ICSD #79500)");
    vproto.push_back("A5B3_mC32_15_e2f_cf;2;16;15;4;14;mC32;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;O5V3;O5V3 (ICSD #15899)");
    vproto.push_back("A5B3_mC32_15_e2f_af;2;16;15;4;14;mC32;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;O5V3;O5V3 (ICSD #32587)");
    vproto.push_back("A8B3_oC22_21_achl_bg;2;11;21;4;8;oC22;a,b/a,c/a,y4,y5,x6,y6,z6;-;O8W3;O8W3 (ICSD #73719)");
    vproto.push_back("AB_oP8_29_a_a;2;8;29;4;9;oP8;a,b/a,c/a,x1,y1,z1,x2,y2,z2;-;O1Pb1;O1Pb1 (ICSD #36250)");
    vproto.push_back("AB_oP4_31_a_a;2;4;31;4;7;oP4;a,b/a,c/a,y1,z1,y2,z2;-;O1Sn1;O1Sn1 (ICSD #20624)");
    vproto.push_back("A5B2_oP14_31_a2b_b;2;14;31;4;14;oP14;a,b/a,c/a,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;D8_{7}/-/-;O5V2 (O5V2)/O5V2/O5V2;O5V2 (ICSD #15984) and O5V2 (shcherbinaite (obsolete), part 3)/O5V2 (ICSD #29140)/O5V2 (ICSD #41030)"); //DX20210428 - added equivalent part 3 prototype (O5V2, shcherbinaite (obsolete), D8_{7}, http://aflow.org/prototype-encyclopedia/A5B2_oP14_31_a2b_b.html)
    vproto.push_back("A2B_oC24_35_abef_de;2;12;35;4;14;oC24;a,b/a,c/a,z1,z2,x3,z3,y4,z4,y5,z5,x6,y6,z6;-;O2Ti1;O2Ti1 (ICSD #97008)");
    vproto.push_back("AB_oC16_36_b_2a;2;8;36;4;10;oC16;a,b/a,c/a,y1,z1,y2,z2,x3,y3,z3;-;O1Sn1;O1Sn1 (ICSD #60619)");
    vproto.push_back("A8B_oC18_38_2cde_a;2;9;38;4;12;oC18;a,b/a,c/a,z1,x2,z2,x3,z3,y4,z4,y5,z5;-;Cu8O1;Cu8O1 (ICSD #62764)");
    vproto.push_back("AB3_oC32_40_c_3c;2;16;40;4;15;oC32;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Cr1O3;Cr1O3 (ICSD #16031)");
    vproto.push_back("A2B3_oF40_43_b_ab;2;10;43;4;10;oF40;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;-;Au2O3 (Ag2O3);Au2O3 (ICSD #8014) and Ag2O3 (part 3)"); //DX20210427 - added equivalent part 3 prototype
    vproto.push_back("AB3_oI8_44_a_ac;2;4;44;4;7;oI8;a,b/a,c/a,z1,z2,x3,z3;-;Li1O3;Li1O3 (ICSD #180565)");
    vproto.push_back("AB_oP2_47_a_e;2;2;47;4;3;oP2;a,b/a,c/a;-;Nb1O1;Nb1O1 (ICSD #95729)");
    vproto.push_back("AB2_oP3_47_a_m;2;3;47;4;4;oP3;a,b/a,c/a,y2;-;Cu1O2;Cu1O2 (ICSD #54126)");
    vproto.push_back("AB4_oP5_47_e_bdr;2;5;47;4;4;oP5;a,b/a,c/a,z4;-;O1Ta4;O1Ta4 (ICSD #76022)");
    vproto.push_back("A5B2_oP7_47_behq_ad;2;14;47;4;4;oP7;a,b/a,c/a,z6;-;O5Ta2;O5Ta2 (ICSD #95462)");
    vproto.push_back("A4B3_oP14_55_gh_bg;2;14;55;4;9;oP14;a,b/a,c/a,x2,y2,x3,y3,x4,y4;-;O4Pb3;O4Pb3 (ICSD #97282)");
    vproto.push_back("AB2_oP12_59_e_2e;2;12;59;4;9;oP12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-/-;Mn1O2/Hf1O2;Mn1O2 (ICSD #54114)/Hf1O2 (ICSD #83863)");
    vproto.push_back("AB2_oP12_61_a_c;2;12;61;4;6;oP12;a,b/a,c/a,x2,y2,z2;-;Hg1O2;Hg1O2 (ICSD #24774)");
    vproto.push_back("AB2_oP12_62_c_2c;2;12;62;4;9;oP12;a,b/a,c/a,x1,z1,x2,z2,x3,z3;-/-/-;Mn1O2/C1Nb2/C1Fe2;Mn1O2 (ICSD #171866)/C1Nb2 (ICSD #31973)/C1Fe2 (ICSD #187138)"); //DX20210120 - added two metal-carbides
    vproto.push_back("A2B3_oC20_63_ac_cf;2;10;63;4;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-/-;Al2O3/Ga2O3;Al2O3 (ICSD #161062)/Ga2O3 (ICSD #162252)");
    vproto.push_back("A3B4_oC28_63_cf_acf;2;14;63;4;9;oC28;a,b/a,c/a,y2,y3,y4,z4,y5,z5;-;Fe3O4;Fe3O4 (ICSD #263010)");
    vproto.push_back("AB2_oC6_65_a_j;2;3;65;4;4;oC6;a,b/a,c/a,y2;-;Ba1O2;Ba1O2 (ICSD #180398)");
    vproto.push_back("AB2_oF12_69_a_h;2;3;69;4;4;oF12;a,b/a,c/a,y2;-;Cu1O2;Cu1O2 (ICSD #85080)");
    vproto.push_back("AB2_oF12_69_a_g;2;3;69;4;4;oF12;a,b/a,c/a,x2;-;Cu1O2;Cu1O2 (ICSD #150886)");
    vproto.push_back("AB_oI8_71_g_e;2;4;71;4;5;oI8;a,b/a,c/a,x1,y2;-/-;O2Rb2 (CsO)/C2Li2;O2Rb2 (ICSD #25528) and CsO (part 3)/C2Li2 (ICSD #25705)"); //DX20210120 - added metal-carbide //DX20210428 - added equivalent part 3 prototype (CsO, http://aflow.org/prototype-encyclopedia/AB_oI8_71_g_i.html)
    vproto.push_back("AB_oI32_72_2j_2j;2;16;72;4;11;oI32;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4;-;Mg1O1;Mg1O1 (ICSD #181462)");
    vproto.push_back("A3B4_oI28_74_ace_hi;2;14;74;4;8;oI28;a,b/a,c/a,z3,y4,z4,x5,z5;-;Fe3O4;Fe3O4 (ICSD #31156)");
    vproto.push_back("AB2_tI24_87_h_2h;2;12;87;4;8;tI24;a,c/a,x1,y1,x2,y2,x3,y3;-;Mn1O2;Mn1O2 (ICSD #20227)");
    vproto.push_back("AB_tI32_88_cd_f;2;16;88;4;5;tI32;a,c/a,x3,y3,z3;-;Ag1O1;Ag1O1 (ICSD #202055)");
    vproto.push_back("A3B_tP16_113_ef_e;2;16;113;4;9;tP16;a,c/a,x1,z1,x2,z2,x3,y3,z3;-;O3W1;O3W1 (ICSD #86144)");
    vproto.push_back("A2B3_tP5_115_g_abc;2;5;115;4;3;tP5;a,c/a,z4;-;Bi2O3;Bi2O3 (ICSD #168808)");
    vproto.push_back("A3B_tP8_127_bg_a;2;8;127;4;3;tP8;a,c/a,x3;-;O3Re1;O3Re1 (ICSD #77680)");
    vproto.push_back("AB_tP4_129_c_a;2;4;129;4;3;tP4;a,c/a,z2;-;Ba1O1;Ba1O1 (ICSD #15301)");
    vproto.push_back("AB2_tP6_129_c_ac;2;6;129;4;4;tP6;a,c/a,z2,z3;-/-;Gd1O2/Eu1O2;Gd1O2 (ICSD #6318)/Eu1O2 (ICSD #631464)");
    vproto.push_back("A3B_tP8_129_cd_c;2;8;129;4;4;tP8;a,c/a,z1,z2;-;O3W1;O3W1 (ICSD #27961)");
    vproto.push_back("A3B_tP16_130_cf_c;2;16;130;4;5;tP16;a,c/a,z1,z2,x3;-;O3W1 (O3W);O3W1 (ICSD #50732) and alpha-WO3 (part 3)"); //DX20210428 - added equivalent part 3 prototype (alpha-WO3, http://aflow.org/prototype-encyclopedia/A3B_tP16_130_cf_c.html)
    vproto.push_back("A2B3_tP10_132_i_be;2;10;132;4;3;tP10;a,c/a,x3;-;Bi2O3;Bi2O3 (ICSD #168807)");
    vproto.push_back("AB_tI4_139_b_a;2;2;139;4;2;tI4;a,c/a;-;O1Pd1;O1Pd1 (ICSD #41617)");
    vproto.push_back("AB3_tI32_140_ab_hk;2;16;140;4;5;tI32;a,c/a,x3,x4,y4;-;K1O3;K1O3 (ICSD #28347)");
    vproto.push_back("A2B_tI24_141_h_c;2;12;141;4;4;tI24;a,c/a,y2,z2;-;O2Ti1;O2Ti1 (ICSD #93098)");
    vproto.push_back("A4B3_tI28_141_cd_ae;2;14;141;4;3;tI28;a,c/a,z4;-;Cu4O3;Cu4O3 (ICSD #77675)");
    vproto.push_back("A6B_hR7_148_f_a;2;7;148;4;5;hR7;a,c/a,x2,y2,z2;-;Hf1O0.1667;Hf1O0.1667 (ICSD #174039)");
    vproto.push_back("A2B3_hP5_150_d_e;2;5;150;4;4;hP5;a,c/a,z1,x2;-;La2O3;La2O3 (ICSD #26864)");
    vproto.push_back("AB_hP6_152_b_a;2;6;152;4;4;hP6;a,c/a,x1,x2;-;Hg1O1;Hg1O1 (ICSD #24062)");
    vproto.push_back("AB3_hR16_155_2c_def;2;16;155;4;9;hR16;a,c/a,x1,x2,y3,y4,x5,y5,z5;-;O1Zr3;O1Zr3 (ICSD #23402)");
    vproto.push_back("AB_hR2_160_a_a;2;2;160;4;4;hR2;a,c/a,x1,x2;-;Bi1O1;Bi1O1 (ICSD #30361)");
    vproto.push_back("A2B3_hR10_160_ab_2b;2;10;160;4;9;hR10;a,c/a,x1,x2,z2,x3,z3,x4,z4;-;Bi2O3;Bi2O3 (ICSD #168810)");
    vproto.push_back("AB6_hP7_162_a_k;2;7;162;4;4;hP7;a,c/a,x2,z2;-;O1Ti6;O1Ti6 (ICSD #20042)");
    vproto.push_back("A3B_hP8_162_k_d;2;8;162;4;4;hP8;a,c/a,x2,z2;-;Ag6O2;Ag6O2 (ICSD #26557)");
    vproto.push_back("AB3_hP8_162_ab_k;2;8;162;4;4;hP8;a,c/a,x3,z3;-;O1Ti3;O1Ti3 (ICSD #23575)");
    vproto.push_back("AB6_hP14_163_c_i;2;14;163;4;5;hP14;a,c/a,x2,y2,z2;-;O1Ti6;O1Ti6 (ICSD #17009)");
    vproto.push_back("AB3_hP16_163_ac_i;2;16;163;4;5;hP16;a,c/a,x3,y3,z3;-;O1Ti3;O1Ti3 (ICSD #20041)");
    vproto.push_back("AB2_hP6_164_d_abd;2;6;164;4;4;hP6;a,c/a,z3,z4;-;Ce1O2;Ce1O2 (ICSD #189287)");
    vproto.push_back("A3B4_hR14_166_acd_ch;2;14;166;4;6;hR14;a,c/a,x2,x3,x5,z5;-;Fe3O4;Fe3O4 (ICSD #92356)");
    vproto.push_back("AB_hR16_166_ab3c_4c;2;16;166;4;9;hR16;a,c/a,x3,x4,x5,x6,x7,x8,x9;-;Mg1O1;Mg1O1 (ICSD #181459)");
    vproto.push_back("AB3_hR8_167_b_e;2;8;167;4;3;hR8;a,c/a,x2;-;O1Zr3;O1Zr3 (ICSD #27023)");
    vproto.push_back("A2B_hP6_186_ab_b;2;6;186;4;5;hP6;a,c/a,z1,z2,z3;-;O2Pt1;O2Pt1 (ICSD #24923)");
    vproto.push_back("A3B_hP12_191_gl_f;2;12;191;4;3;hP12;a,c/a,x3;-;O3W1 (O3W);O3W1 (ICSD #32001) and Hexagonal WO3 (part 3)"); //DX20210428 - added equivalent part 3 prototype (hexagonal WO3, http://aflow.org/prototype-encyclopedia/A3B_hP12_191_gl_f.html)
    vproto.push_back("AB_hP8_194_f_f;2;8;194;4;4;hP8;a,c/a,z1,z2;-/-;Li2O2/Li2O2;Li2O2 (ICSD #24143)/Li2O2 (ICSD #50658)");
    vproto.push_back("A3B4_cF56_216_abe_2e;2;14;216;4;4;cF56;a,x3,x4,x5;-;Fe3O4;Fe3O4 (ICSD #65338)");
    vproto.push_back("AB2_cI24_217_c_abc;2;12;217;4;3;cI24;a,x3,x4;-;O1Ta2;O1Ta2 (ICSD #28387)");
    vproto.push_back("AB_cP12_223_c_d;2;12;223;4;1;cP12;a;-;Mg1O1;Mg1O1 (ICSD #181465)");
    vproto.push_back("A4B3_cP14_223_e_c;2;14;223;4;1;cP14;a;-;O4Pt3;O4Pt3 (ICSD #30444)");
    vproto.push_back("A2B3_cP10_224_b_d;2;10;224;4;1;cP10;a;D5_{5};Ag2O3 (Mg3P2);Ag2O3 (ICSD #15999) and Mg3P2 (part 3)"); //DX20210428 - add equivalent part 3 prototype (Mg3P2, D5_{5}, http://aflow.org/prototype-encyclopedia/A3B2_cP10_224_d_b.html, part 3)
    // -------------------------------------------------------------------------
    // ternaries
    // -------------------------------------------------------------------------
    vproto.push_back("AB3C_aP5_1_a_3a_a;3;5;1;5;21;aP5;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Cu1O3V1;Cu1O3V1 (ICSD #9414)");
    vproto.push_back("ABC3_aP10_1_2a_2a_6a;3;10;1;5;36;aP10;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Ca1Fe1O3;Ca1Fe1O3 (ICSD #92339)");
    vproto.push_back("A2B4C_aP14_1_4a_8a_2a;3;14;1;5;48;aP14;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Al2O4Pb1;Al2O4Pb1 (ICSD #33532)");
    vproto.push_back("A2BC4_aP14_1_4a_2a_8a;3;14;1;5;48;aP14;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Cr2Ni1O4;Cr2Ni1O4 (ICSD #280061)");
    vproto.push_back("AB4C3_aP16_1_2a_8a_6a;3;16;1;5;54;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;Co1Na4O3;Co1Na4O3 (ICSD #10473)");
    vproto.push_back("AB6C2_aP9_2_a_3i_i;3;9;2;5;18;aP9;a,b/a,c/a,alpha,beta,gamma,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Cu1O6V2;Cu1O6V2 (ICSD #28151)");
    vproto.push_back("AB3C_aP10_2_i_3i_i;3;10;2;5;21;aP10;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-/-;Mn1O3Sn1/Hg1O3V1;Mn1O3Sn1 (ICSD #29203)/Hg1O3V1 (ICSD #82242)");
    vproto.push_back("AB3C6_aP10_2_a_bi_3i;3;10;2;5;18;aP10;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Au1Gd3O6;Au1Gd3O6 (ICSD #411500)");
    vproto.push_back("AB4C_aP12_2_i_4i_i;3;12;2;5;24;aP12;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Cu1O4W1;Cu1O4W1 (ICSD #4189)");
    vproto.push_back("A3B8C2_aP13_2_ai_4i_i;3;13;2;5;24;aP13;a,b/a,c/a,alpha,beta,gamma,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Cu3O8V2;Cu3O8V2 (ICSD #27184)");
    vproto.push_back("AB3C3_aP14_2_i_3i_3i;3;14;2;5;27;aP14;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Co2Na6O6;Co2Na6O6 (ICSD #99580)");
    vproto.push_back("A4BC2_aP14_2_4i_i_abcf;3;14;2;5;21;aP14;a,b/a,c/a,alpha,beta,gamma,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;O4Pb1Pt2;O4Pb1Pt2 (ICSD #59657)");
    vproto.push_back("A4BC3_aP16_2_4i_i_3i;3;16;2;5;30;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-/-;Cs4Fe1O3/C4K1N3;Cs4Fe1O3 (ICSD #423336)/C4K1N3 (ICSD #77046)"); //DX20210221 - added carbo-nitride
    vproto.push_back("ABC2_mP16_4_2a_2a_4a;3;16;4;5;28;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Ca1Cu1O2;Ca1Cu1O2 (ICSD #84868)");
    vproto.push_back("AB6C2_mC18_5_a_3c_c;3;9;5;5;17;mC18;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Cu1O6V2;Cu1O6V2 (ICSD #21067)");
    vproto.push_back("A6B2C_mC18_5_3c_c_a;3;9;5;5;17;mC18;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;O6V2Zn1;O6V2Zn1 (ICSD #26998)");
    vproto.push_back("ABC4_mC24_5_c_c_4c;3;12;5;5;22;mC24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Fe1Nb1O4;Fe1Nb1O4 (ICSD #429)");
    vproto.push_back("A8B3C2_mC26_5_4c_ac_c;3;13;5;5;23;mC26;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;O8Pb3V2;O8Pb3V2 (ICSD #29359)");
    vproto.push_back("AB4C2_mC28_5_c_4c_abc;3;14;5;5;24;mC28;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Mo1O4Tl2;Mo1O4Tl2 (ICSD #280608)");
    vproto.push_back("AB3C_mP5_6_b_2ab_a;3;5;6;5;14;mP5;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;Ba1O3Ti1;Ba1O3Ti1 (ICSD #186460)");
    vproto.push_back("A2BC3_mC12_8_2a_a_3a;3;6;8;5;16;mC12;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Ca2Co1O3;Ca2Co1O3 (ICSD #95439)");
    vproto.push_back("AB5C5_mC22_8_a_5a_5a;3;11;8;5;26;mC22;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11;-;Bi1Li5O5;Bi1Li5O5 (ICSD #203031)");
    vproto.push_back("ABC2_mC24_8_ab_ab_2a2b;3;12;8;5;24;mC24;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Li1Ni1O2;Li1Ni1O2 (ICSD #164853)");
    vproto.push_back("A2B6C5_mC26_8_2a_2a2b_a2b;3;13;8;5;26;mC26;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Fe2K6O5;Fe2K6O5 (ICSD #174312)");
    vproto.push_back("A6B2C5_mC26_8_2a2b_2a_a2b;3;13;8;5;26;mC26;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Cs6Fe2O5;Cs6Fe2O5 (ICSD #174313)");
    vproto.push_back("ABC3_mC20_9_a_a_3a;3;10;9;5;19;mC20;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Bi1Fe1O3;Bi1Fe1O3 (ICSD #247765)");
    vproto.push_back("A2B3C_mC24_9_2a_3a_a;3;12;9;5;22;mC24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Li2O3Zr1;Li2O3Zr1 (ICSD #31941)");
    vproto.push_back("AB4C3_mC32_9_a_4a_3a;3;16;9;5;28;mC32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-/-;Fe1Na4O3/Co1Na4O3;Fe1Na4O3 (ICSD #1410)/Co1Na4O3 (ICSD #10501)");
    vproto.push_back("A4B3C_mC32_9_4a_3a_a;3;16;9;5;28;mC32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Na4O3Sn1;Na4O3Sn1 (ICSD #49624)");
    vproto.push_back("ABC3_mP10_11_a_e_d2e;3;10;11;5;10;mP10;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5;-;Bi1Fe1O3;Bi1Fe1O3 (ICSD #162834)");
    vproto.push_back("AB4C_mP12_11_e_aef_e;3;12;11;5;13;mP12;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;Ba1O4W1;Ba1O4W1 (ICSD #155516)");
    vproto.push_back("A2BC4_mP14_11_2e_e_4e;3;14;11;5;18;mP14;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Al2Ca1O4;Al2Ca1O4 (ICSD #172780)");
    vproto.push_back("ABC2_mC8_12_a_c_i;3;4;12;5;6;mC8;a,b/a,c/a,beta,x3,z3;-/-/-;Cu1Na1O2/Mn1Na1O2/Cu1Li1O2;Cu1Na1O2 (ICSD #15098)/Mn1Na1O2 (ICSD #16270)/Cu1Li1O2 (ICSD #74978)");
    vproto.push_back("ABC2_mC8_12_c_a_i;3;4;12;5;6;mC8;a,b/a,c/a,beta,x3,z3;-/-;Na1Ni1O2/Cu1Mn1O2;Na1Ni1O2 (ICSD #26609)/Cu1Mn1O2 (ICSD #30379)");
    vproto.push_back("ABC2_mC8_12_a_d_i;3;4;12;5;6;mC8;a,b/a,c/a,beta,x3,z3;-/-;Cd1Hg1O2/Ag1Cu1O2;Cd1Hg1O2 (ICSD #74848)/Ag1Cu1O2 (ICSD #95089)");
    vproto.push_back("AB2C_mC8_12_a_i_b;3;4;12;5;6;mC8;a,b/a,c/a,beta,x3,z3;-;Na1O2V1;Na1O2V1 (ICSD #420136)");
    vproto.push_back("AB2C_mC8_12_a_i_c;3;4;12;5;6;mC8;a,b/a,c/a,beta,x3,z3;-;Na1O2V1;Na1O2V1 (ICSD #420138)");
    vproto.push_back("AB2C2_mC10_12_a_i_i;3;5;12;5;8;mC10;a,b/a,c/a,beta,x2,z2,x3,z3;-;Cu1Li2O2/C1K2N2/C1N2Na2;Cu1Li2O2 (ICSD #174134)/C1K2N2 (ICSD #411094)/C1N2Na2 (ICSD #411341)"); //DX20210221 - added two carbo-nitrides
    vproto.push_back("ABC2_mC16_12_ad_i_2i;3;8;12;5;10;mC16;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5;-;Li1Mo1O2;Li1Mo1O2 (ICSD #165326)");
    vproto.push_back("AB4C4_mC18_12_a_2i_j;3;9;12;5;11;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,y4,z4;-;Ir1K4O4;Ir1K4O4 (ICSD #47223)");
    vproto.push_back("A2B3C4_mC18_12_i_ai_2i;3;9;12;5;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-;Cu2Li3O4;Cu2Li3O4 (ICSD #66509)");
    vproto.push_back("A2B5C2_mC18_12_i_a2i_i;3;9;12;5;12;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5;-;K2O5Ti2 (K2O5Ti2);K2O5Ti2 (ICSD #36097) and K2O5Ti2 (part 3)"); //DX20210427 - added equivalent part 3 prototype
    vproto.push_back("A4BC4_mC18_12_2i_a_j;3;9;12;5;11;mC18;a,b/a,c/a,beta,x2,z2,x3,z3,x4,y4,z4;-;Cs4Ir1O4;Cs4Ir1O4 (ICSD #72287)");
    vproto.push_back("AB3C_mC20_12_i_ghi_e;3;10;12;5;10;mC20;a,b/a,c/a,beta,y2,y3,x4,z4,x5,z5;-;Ba1O3Pb1;Ba1O3Pb1 (ICSD #94313)");
    vproto.push_back("AB8C2_mC22_12_b_ac3i_i;3;11;12;5;12;mC22;a,b/a,c/a,beta,x4,z4,x5,z5,x6,z6,x7,z7;-;Mo1O8V2;Mo1O8V2 (ICSD #28471)");
    vproto.push_back("A2B8C_mC22_12_i_2ij_a;3;11;12;5;13;mC22;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;Mo2O8Zr1;Mo2O8Zr1 (ICSD #280433)");
    vproto.push_back("A7B2C2_mC22_12_aij_i_h;3;11;12;5;12;mC22;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,y5,z5;-;O7V2Zn2;O7V2Zn2 (ICSD #250002)");
    vproto.push_back("A3BC8_mC24_12_ci_a_2ij;3;12;12;5;13;mC24;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,y6,z6;-/-/-/-;Cr3K1O8/Cr3Na1O8/Cr3In1O8/Cr3K1O8;Cr3K1O8 (ICSD #15898)/Cr3Na1O8 (ICSD #82620)/Cr3In1O8 (ICSD #155504)/Cr3K1O8 (ICSD #246524)");
    vproto.push_back("ABC4_mC24_12_i_i_4i;3;12;12;5;16;mC24;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Fe1Nb1O4 (AlNbO4);Fe1Nb1O4 (ICSD #14016) and AlNbO4 (part 3)"); //DX20210428 - added equivalent part 3 prototype (AlNbO4, http://aflow.org/prototype-encyclopedia/ABC4_mC24_12_i_i_4i.html)
    vproto.push_back("AB4C_mC24_12_i_2ij_g;3;12;12;5;14;mC24;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;Al1O4W1;Al1O4W1 (ICSD #4164)");
    vproto.push_back("A2BC3_mC24_12_acg_h_ij;3;12;12;5;11;mC24;a,b/a,c/a,beta,y3,y4,x5,z5,x6,y6,z6;-;Li2Ni1O3;Li2Ni1O3 (ICSD #153094)");
    vproto.push_back("A4B7C_mC24_12_2i_b3i_a;3;12;12;5;14;mC24;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;La4O7Pd1;La4O7Pd1 (ICSD #65032)");
    vproto.push_back("A5B6C_mC24_12_dgh_ij_a;3;12;12;5;11;mC24;a,b/a,c/a,beta,y3,y4,x5,z5,x6,y6,z6;-/-;Li5O6Re1/Na5O6Re1;Li5O6Re1 (ICSD #38381)/Na5O6Re1 (ICSD #38382)");
    vproto.push_back("A3B8C_mC24_12_ci_2ij_a;3;12;12;5;13;mC24;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;Cr3O8Tl1;Cr3O8Tl1 (ICSD #155505)");
    vproto.push_back("AB3C8_mC24_12_a_ci_2ij;3;12;12;5;13;mC24;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;Ag1Cr3O8;Ag1Cr3O8 (ICSD #155507)");
    vproto.push_back("A3BC8_mC24_12_di_a_2ij;3;12;12;5;13;mC24;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;Cr3Li1O8;Cr3Li1O8 (ICSD #155851)");
    vproto.push_back("A4BC8_mC26_12_2i_a_4i;3;13;12;5;16;mC26;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-/-;Cr8K2O16/Ir1K0.25O2;Cr8K2O16 (ICSD #491)/Ir1K0.25O2 (ICSD #80336)");
    vproto.push_back("A2B3C8_mC26_12_i_ah_2ij;3;13;12;5;14;mC26;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-/-;Cu2Mn3O8/Cd2Mn3O8;Cu2Mn3O8 (ICSD #971)/Cd2Mn3O8 (ICSD #16957)");
    vproto.push_back("A3B8C2_mC26_12_ai_2ij_i;3;13;12;5;15;mC26;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;Ca3O8V2;Ca3O8V2 (ICSD #412273)");
    vproto.push_back("A6B2C5_mC26_12_ghi_i_aj;3;13;12;5;13;mC26;a,b/a,c/a,beta,y2,y3,x4,z4,x5,z5,x6,y6,z6;-;Cs6Fe2O5;Cs6Fe2O5 (ICSD #73134)");
    vproto.push_back("A4B3C6_mC26_12_2i_ai_hj;3;13;12;5;14;mC26;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;K4Ni3O6;K4Ni3O6 (ICSD #426450)");
    vproto.push_back("A4B7C2_mC26_12_aci_bij_i;3;13;12;5;13;mC26;a,b/a,c/a,beta,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;K4O7V2;K4O7V2 (ICSD #250388)");
    vproto.push_back("AB4C8_mC26_12_a_2i_4i;3;13;12;5;16;mC26;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Ba2Mn8O16;Ba2Mn8O16 (ICSD #62096)");
    vproto.push_back("A4B8C_mC26_12_2i_4i_a;3;13;12;5;16;mC26;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Ga4O8Ti1;Ga4O8Ti1 (ICSD #155638)");
    vproto.push_back("AB3C3_mC28_12_i_ghi_ij;3;14;12;5;15;mC28;a,b/a,c/a,beta,y1,y2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;Ga2K6O6;Ga2K6O6 (ICSD #2269)");
    vproto.push_back("AB2C4_mC28_12_i_2i_4i;3;14;12;5;18;mC28;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Hg1Ho2O4;Hg1Ho2O4 (ICSD #69734)");
    vproto.push_back("A2BC4_mC28_12_2i_i_2ij;3;14;12;5;17;mC28;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;K2Mo1O4;K2Mo1O4 (ICSD #16154)");
    vproto.push_back("A2B4C_mC28_12_2i_4i_i;3;14;12;5;18;mC28;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;Bi2O4Sr1;Bi2O4Sr1 (ICSD #80668)");
    vproto.push_back("A3B10C2_mC30_12_ah_3ij_i;3;15;12;5;16;mC30;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;La3O10Os2;La3O10Os2 (ICSD #10104)");
    vproto.push_back("A10B3C2_mC30_12_3ij_ah_i;3;15;12;5;16;mC30;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;O10Ru3Sr2;O10Ru3Sr2 (ICSD #50707)");
    vproto.push_back("AB5C2_mC32_12_i_g2ij_2i;3;16;12;5;18;mC32;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-/-;Cr1O5Pb2/Cr1O5Pb2;Cr1O5Pb2 (ICSD #29269)/Cr1O5Pb2 (ICSD #201562)");
    vproto.push_back("A5B2C_mC32_12_g2ij_2i_i;3;16;12;5;18;mC32;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;-;O5Pb2W1;O5Pb2W1 (ICSD #61399)");
    vproto.push_back("AB4C_mP12_13_e_2g_f;3;12;13;5;12;mP12;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-/-;Mo1O4Zn1/Nd1O4Ta1;Mo1O4Zn1 (ICSD #236416)/Nd1O4Ta1 (ICSD #415430)");
    vproto.push_back("AB3C4_mP16_13_e_e2f_2g;3;16;13;5;14;mP16;a,b/a,c/a,beta,y1,y2,y3,y4,x5,y5,z5,x6,y6,z6;-;Bi1Na3O4;Bi1Na3O4 (ICSD #10319)");
    vproto.push_back("A2B2C3_mP14_14_e_e_ae;3;14;14;5;13;mP14;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Cd2K2O3;Cd2K2O3 (ICSD #16223)");
    vproto.push_back("AB2C_mP16_14_e_2e_e;3;16;14;5;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Li1O2Y1;Li1O2Y1 (ICSD #45511)");
    vproto.push_back("ABC2_mP16_14_e_e_2e;3;16;14;5;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-/E0_{6};Ce1Li1O2/C1Hg1N2/KNO2/MnHO2;Ce1Li1O2 (ICSD #47116)/C1Hg1N2 (ICSD #412278)/KNO2 III/Manganite (gamma-MnO(OH), E06)"); //DX20210221 - added carbo-nitride //DX20210427 - added KNO2 III and Manganite (gamma-MnO(OH), E06) from part 3
    vproto.push_back("ABC2_mC16_15_e_e_f;3;8;15;5;9;mC16;a,b/a,c/a,beta,y1,y2,x3,y3,z3;-/-/-/F5_{a};Bi1Na1O2/Bi1K1O2/Bi1Cs1O2/FeKS2;Bi1Na1O2 (ICSD #10317)/Bi1K1O2 (ICSD #16239)/Bi1Cs1O2 (ICSD #406564)/KFeS2 (F5a)"); //DX20210427 - added KFeS2 (F5a) from part 3
    vproto.push_back("AB2C_mC16_15_e_f_e;3;8;15;5;9;mC16;a,b/a,c/a,beta,y1,y2,x3,y3,z3;-;Bi1O2Rb1;Bi1O2Rb1 (ICSD #407208)");
    vproto.push_back("ABC3_mC20_15_a_e_ef;3;10;15;5;9;mC20;a,b/a,c/a,beta,y2,y3,x4,y4,z4;-;Al1La1O3;Al1La1O3 (ICSD #180415)");
    vproto.push_back("A2B2C_mC20_15_ac_f_e;3;10;15;5;8;mC20;a,b/a,c/a,beta,y3,x4,y4,z4;-;Cu2O2Pb1;Cu2O2Pb1 (ICSD #400657)");
    vproto.push_back("ABC4_mC24_15_a_e_2f;3;12;15;5;11;mC24;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;-;Hg1Mo1O4;Hg1Mo1O4 (ICSD #2533)");
    vproto.push_back("ABC4_mC24_15_e_e_2f;3;12;15;5;12;mC24;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;La1Nb1O4;La1Nb1O4 (ICSD #10123)");
    vproto.push_back("AB4C_mC24_15_e_2f_e;3;12;15;5;12;mC24;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4;-;Nb1O4Y1;Nb1O4Y1 (ICSD #20335)");
    vproto.push_back("AB4C_mC24_15_a_2f_e;3;12;15;5;11;mC24;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;-;Hg1O4W1;Hg1O4W1 (ICSD #169668)");
    vproto.push_back("A2BC4_mC28_15_f_e_2f;3;14;15;5;14;mC28;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Hg2Mo1O4;Hg2Mo1O4 (ICSD #90084)");
    vproto.push_back("A3B4C_mC32_15_cf_2f_e;3;16;15;5;14;mC32;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Ag3O4V1;Ag3O4V1 (ICSD #249417)");
    vproto.push_back("A2B5C_mC32_15_f_e2f_c;3;16;15;5;14;mC32;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Fe2O5Ti1;Fe2O5Ti1 (ICSD #24416)");
    vproto.push_back("ABC2_mC32_15_2e_2e_2f;3;16;15;5;14;mC32;a,b/a,c/a,beta,y1,y2,y3,y4,x5,y5,z5,x6,y6,z6;-;Er1Na1O2;Er1Na1O2 (ICSD #2739)");
    vproto.push_back("A3B4C_oP16_31_ab_2ab_a;3;16;31;5;17;oP16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;-;Li3O4V1;Li3O4V1 (ICSD #19002)");
    vproto.push_back("ABC3_oC20_36_a_a_3a;3;10;36;5;13;oC20;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5;-;Ca1Ir1O3;Ca1Ir1O3 (ICSD #159030)");
    vproto.push_back("AB3C_oC20_36_a_ab_a;3;10;36;5;12;oC20;a,b/a,c/a,y1,z1,y2,z2,y3,z3,x4,y4,z4;-;Mg1O3V1;Mg1O3V1 (ICSD #284)");
    vproto.push_back("A2BC2_oC20_36_2a_a_2a;3;10;36;5;13;oC20;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5;-;Na2Ni1O2;Na2Ni1O2 (ICSD #14159)");
    vproto.push_back("AB4C_oC24_36_a_4a_a;3;12;36;5;15;oC24;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6;-;La1O4Ta1;La1O4Ta1 (ICSD #97688)");
    vproto.push_back("AB5C2_oC32_36_a_5a_2a;3;16;36;5;19;oC32;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8;-;Li1O5V2;Li1O5V2 (ICSD #23817)");
    vproto.push_back("ABC3_oC10_38_a_b_ae;3;5;38;5;8;oC10;a,b/a,c/a,z1,z2,z3,y4,z4;-;K1Nb1O3;K1Nb1O3 (ICSD #9533)");
    vproto.push_back("AB3C_oC10_38_b_ad_a;3;5;38;5;8;oC10;a,b/a,c/a,z1,z2,z3,y4,z4;-/-;Ba1O3Ti1/Ba1O3Ti1;Ba1O3Ti1 (ICSD #154346)/Ba1O3Ti1 (ICSD #161341)");
    vproto.push_back("AB3C_oC10_38_a_ae_b;3;5;38;5;8;oC10;a,b/a,c/a,z1,z2,z3,y4,z4;-;Ba1O3Ti1;Ba1O3Ti1 (ICSD #161419)");
    vproto.push_back("ABC_oC12_40_b_a_b;3;6;40;5;8;oC12;a,b/a,c/a,z1,y2,z2,y3,z3;-;Cs1Cu1O1;Cs1Cu1O1 (ICSD #37079)");
    vproto.push_back("A2B4C_oC28_40_c_abc_b;3;14;40;5;14;oC28;a,b/a,c/a,z1,y2,z2,y3,z3,x4,y4,z4,x5,y5,z5;-;Al2O4Pb1;Al2O4Pb1 (ICSD #80128)");
    vproto.push_back("A3BC_oI20_46_a2b_a_b;3;10;46;5;11;oI20;a,b/a,c/a,z1,z2,y3,z3,y4,z4,y5,z5;-;O3Sr1Ti1;O3Sr1Ti1 (ICSD #182248)");
    vproto.push_back("ABC2_oP4_47_a_g_s;3;4;47;5;4;oP4;a,b/a,c/a,z3;-;Au1K1O2;Au1K1O2 (ICSD #15115)");
    vproto.push_back("A3BC_oP5_47_bce_h_a;3;5;47;5;3;oP5;a,b/a,c/a;-;O3Pb1Ti1;O3Pb1Ti1 (ICSD #27949)");
    vproto.push_back("AB4C_oP12_50_a_m_c;3;12;50;5;6;oP12;a,b/a,c/a,x3,y3,z3;-;Ce1O4Zr1;Ce1O4Zr1 (ICSD #164740)");
    vproto.push_back("A4BC2_oP14_55_gh_a_h;3;14;55;5;9;oP14;a,b/a,c/a,x2,y2,x3,y3,x4,y4;-;O4Pb1Sr2;O4Pb1Sr2 (ICSD #4418)");
    vproto.push_back("A2B4C_oP14_55_h_gh_a;3;14;55;5;9;oP14;a,b/a,c/a,x2,y2,x3,y3,x4,y4;-;Cd2O4Sn1;Cd2O4Sn1 (ICSD #9010)");
    vproto.push_back("AB2C_oP8_58_c_g_a;3;8;58;5;5;oP8;a,b/a,c/a,x3,y3;-;Li1O2Ru1;Li1O2Ru1 (ICSD #48007)");
    vproto.push_back("A2BC3_oP12_59_e_a_be;3;12;59;5;9;oP12;a,b/a,c/a,z1,z2,y3,z3,y4,z4;-;Cu2Mg1O3;Cu2Mg1O3 (ICSD #4202)");
    vproto.push_back("A2B3C_oP12_59_e_be_a;3;12;59;5;9;oP12;a,b/a,c/a,z1,z2,y3,z3,y4,z4;-;Cu2O3Sr1;Cu2O3Sr1 (ICSD #95957)");
    vproto.push_back("AB5C2_oP16_59_b_a2e_e;3;16;59;5;11;oP16;a,b/a,c/a,z1,z2,y3,z3,y4,z4,y5,z5;-;Na1O5V2;Na1O5V2 (ICSD #59345)");
    vproto.push_back("AB2C_oP16_59_e_2e_e;3;16;59;5;11;oP16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-;Li1O2Tb1;Li1O2Tb1 (ICSD #21013)");
    vproto.push_back("AB5C2_oP16_59_a_b2e_e;3;16;59;5;11;oP16;a,b/a,c/a,z1,z2,y3,z3,y4,z4,y5,z5;-;Li1O5V2;Li1O5V2 (ICSD #88643)");
    vproto.push_back("ABC2_oP16_59_e_e_2e;3;16;59;5;11;oP16;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4;-;Gd1Li1O2;Gd1Li1O2 (ICSD #27769)");
    vproto.push_back("ABC2_oP16_62_c_c_2c;3;16;62;5;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;-/E0_{2};Eu1Li1O2/AlHO2;Eu1Li1O2 (ICSD #422560)/Diaspore (AlOOH, E02)"); //DX20210427 - added Diaspore (AlOOH, E02) from part 3
    vproto.push_back("ABC_oC12_63_c_a_c;3;6;63;5;5;oC12;a,b/a,c/a,y2,y3;-;Cs1Cu1O1;Cs1Cu1O1 (ICSD #40161)");
    vproto.push_back("ABC2_oC16_63_c_c_f;3;8;63;5;7;oC16;a,b/a,c/a,y1,y2,y3,z3;-/-/-;Cu1K1O2/Cs1Cu1O2/Au1Na1O2;Cu1K1O2 (ICSD #15095)/Cs1Cu1O2 (ICSD #15097)/Au1Na1O2 (ICSD #409547)");
    vproto.push_back("ABC2_oC16_63_c_a_g;3;8;63;5;6;oC16;a,b/a,c/a,y2,x3,y3;-;Ba1Ni1O2;Ba1Ni1O2 (ICSD #15760)");
    vproto.push_back("ABC3_oC20_63_a_c_cf;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;Ir1Na1O3;Ir1Na1O3 (ICSD #261371)");
    vproto.push_back("AB4C_oC24_63_c_acf_c;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,y4,y5,z5;-;Bi1O4Re1;Bi1O4Re1 (ICSD #10481)");
    vproto.push_back("A2B3C_oC24_63_g_ce_c;3;12;63;5;8;oC24;a,b/a,c/a,y1,y2,x3,x4,y4;-;Cs2O3Pb1;Cs2O3Pb1 (ICSD #62140)");
    vproto.push_back("ABC4_oC24_63_c_a_fg;3;12;63;5;8;oC24;a,b/a,c/a,y2,y3,z3,x4,y4;-/-;Cr1Cu1O4/Cr1Hg1O4;Cr1Cu1O4 (ICSD #60825)/Cr1Hg1O4 (ICSD #416147)");
    vproto.push_back("AB4C2_oC28_63_c_acf_f;3;14;63;5;9;oC28;a,b/a,c/a,y2,y3,y4,z4,y5,z5;-;Ca1O4Tl2;Ca1O4Tl2 (ICSD #79371)");
    vproto.push_back("A2BC4_oC28_63_f_c_acf;3;14;63;5;9;oC28;a,b/a,c/a,y2,y3,y4,z4,y5,z5;-;Al2Mg1O4;Al2Mg1O4 (ICSD #161057)");
    vproto.push_back("ABC_oC24_64_d_f_f;3;12;64;5;8;oC24;a,b/a,c/a,x1,y2,z2,y3,z3;-;K1Li1O1;K1Li1O1 (ICSD #30964)");
    vproto.push_back("AB2C4_oC28_64_a_f_ef;3;14;64;5;8;oC28;a,b/a,c/a,y2,y3,z3,y4,z4;-/-;Cu1La2O4/Cu1Nd2O4;Cu1La2O4 (ICSD #67836)/Cu1Nd2O4 (ICSD #86754)");
    vproto.push_back("AB4C2_oC28_64_a_ef_f;3;14;64;5;8;oC28;a,b/a,c/a,y2,y3,z3,y4,z4;-;Cu1O4Pr2;Cu1O4Pr2 (ICSD #91073)");
    vproto.push_back("ABC3_oC10_65_a_c_bf;3;5;65;5;3;oC10;a,b/a,c/a;-;Na1Nb1O3;Na1Nb1O3 (ICSD #28564)");
    vproto.push_back("AB4C_oC12_65_a_hi_c;3;6;65;5;5;oC12;a,b/a,c/a,x3,y4;-;Co1O4Re1;Co1O4Re1 (ICSD #72872)");
    vproto.push_back("A2B3C_oC12_65_i_ai_c;3;6;65;5;5;oC12;a,b/a,c/a,y3,y4;-;Cu2O3Sr1;Cu2O3Sr1 (ICSD #99041)");
    vproto.push_back("A2B3C_oC12_65_j_bj_a;3;6;65;5;5;oC12;a,b/a,c/a,y3,y4;-;Li2O3Pr1 (Li2O3Pr);Li2O3Pr1 (ICSD #154704) and Li2O3Pr (part 3)"); //DX20210427 - added equivalent part 3 proto
    vproto.push_back("AB3C4_oC16_65_b_ai_p;3;8;65;5;6;oC16;a,b/a,c/a,y3,x4,y4;-;Ba1Cu3O4;Ba1Cu3O4 (ICSD #65881)");
    vproto.push_back("AB2C_oC16_65_i_n_ad;3;8;65;5;6;oC16;a,b/a,c/a,y3,y4,z4;-;O1Ti2Zr1;O1Ti2Zr1 (ICSD #9389)");
    vproto.push_back("AB3C4_oC16_65_c_ai_p;3;8;65;5;6;oC16;a,b/a,c/a,y3,x4,y4;-;Ba1Cu3O4;Ba1Cu3O4 (ICSD #83079)");
    vproto.push_back("AB3C4_oC16_65_d_ai_p;3;8;65;5;6;oC16;a,b/a,c/a,y3,x4,y4;-;Ba1Cu3O4;Ba1Cu3O4 (ICSD #89232)");
    vproto.push_back("AB6C3_oC20_65_b_jp_af;3;10;65;5;6;oC20;a,b/a,c/a,y4,x5,y5;-;Mn1O6Pt3;Mn1O6Pt3 (ICSD #35337)");
    vproto.push_back("A3B5C2_oC20_65_ai_b2i_j;3;10;65;5;7;oC20;a,b/a,c/a,y3,y4,y5,y6;-;Cu6O10Sr4;Cu6O10Sr4 (ICSD #50089)");
    vproto.push_back("AB3C_oC20_65_i_acm_j;3;10;65;5;6;oC20;a,b/a,c/a,y3,y4,z5;-;Mg1O3V1;Mg1O3V1 (ICSD #15927)");
    vproto.push_back("A3B7C_oC22_65_af_cjp_b;3;11;65;5;6;oC22;a,b/a,c/a,y5,x6,y6;-;La3O7Ta1;La3O7Ta1 (ICSD #168909)");
    vproto.push_back("A2BC4_oC28_65_ij_ac_ijo;3;14;65;5;9;oC28;a,b/a,c/a,y3,y4,y5,y6,x7,z7;-;La2Ni1O4;La2Ni1O4 (ICSD #78338)");
    vproto.push_back("A2BC_oC16_67_m_g_a;3;8;67;5;6;oC16;a,b/a,c/a,z2,y3,z3;-;Ba2Na1O1;Ba2Na1O1 (ICSD #411905)");
    vproto.push_back("A2B3C2_oF28_69_h_ah_g;3;7;69;5;6;oF28;a,b/a,c/a,x2,y3,y4;-;Cu2O3Sr2;Cu2O3Sr2 (ICSD #68676)");
    vproto.push_back("AB4C2_oF28_69_a_cg_g;3;7;69;5;5;oF28;a,b/a,c/a,x3,x4;-;Ni1O4Pr2;Ni1O4Pr2 (ICSD #81577)");
    vproto.push_back("A4BC_oF48_70_h_b_a;3;12;70;5;6;oF48;a,b/a,c/a,x3,y3,z3;-;O4Tb1V1;O4Tb1V1 (ICSD #88369)");
    vproto.push_back("A2BC4_oF56_70_d_a_h;3;14;70;5;6;oF56;a,b/a,c/a,x3,y3,z3;-;Cr2Mg1O4;Cr2Mg1O4 (ICSD #290589)");
    vproto.push_back("AB2C4_oF56_70_a_g_h;3;14;70;5;7;oF56;a,b/a,c/a,z2,x3,y3,z3;-;Mo1Na2O4;Mo1Na2O4 (ICSD #151971)");
    vproto.push_back("A2B2C_oI10_71_e_f_a;3;5;71;5;5;oI10;a,b/a,c/a,x2,x3;-/-/-/-;Li2O2Pd1/Na2O2Pt1/Ag2O2Pd1/Li2O2Pd1;Li2O2Pd1 (ICSD #19007)/Na2O2Pt1 (ICSD #25018)/Ag2O2Pd1 (ICSD #51498)/Li2O2Pd1 (ICSD #61199)");
    vproto.push_back("AB2C2_oI10_71_a_e_f;3;5;71;5;5;oI10;a,b/a,c/a,x2,x3;-;Cu1Li2O2;Cu1Li2O2 (ICSD #25001)");
    vproto.push_back("AB2C3_oI12_71_b_e_af;3;6;71;5;5;oI12;a,b/a,c/a,x3,x4;-;Ba1Mn2O3;Ba1Mn2O3 (ICSD #10038)");
    vproto.push_back("ABC2_oI16_71_e_e_2f;3;8;71;5;7;oI16;a,b/a,c/a,x1,x2,x3,x4;-;Li1Mn1O2;Li1Mn1O2 (ICSD #15768)");
    vproto.push_back("A2B4C3_oI18_71_f_n_ae;3;9;71;5;7;oI18;a,b/a,c/a,x2,x3,x4,y4;-;Na2O4Pd3;Na2O4Pd3 (ICSD #6157)");
    vproto.push_back("A3B5C2_oI20_71_ae_def_e;3;10;71;5;7;oI20;a,b/a,c/a,x3,x4,x5,x6;-;Cu3O5Sr2;Cu3O5Sr2 (ICSD #416904)");
    vproto.push_back("A5B2C7_oI28_71_aef_g_dhn;3;14;71;5;9;oI28;a,b/a,c/a,x3,x4,y5,y6,x7,y7;-;Au5O2Rb7;Au5O2Rb7 (ICSD #91309)");
    vproto.push_back("AB3C2_oI24_72_c_bf_j;3;12;72;5;6;oI24;a,b/a,c/a,x3,x4,y4;-;Ag1Na3O2;Ag1Na3O2 (ICSD #16919)");
    vproto.push_back("A3BC2_oI24_72_ce_b_j;3;12;72;5;5;oI24;a,b/a,c/a,x4,y4;-;Ag3Li1O2;Ag3Li1O2 (ICSD #4204)");
    vproto.push_back("ABC2_oI32_72_j_f_2j;3;16;72;5;10;oI32;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;-;Bi1Li1O2;Bi1Li1O2 (ICSD #25385)");
    vproto.push_back("A2BC_oI16_74_f_e_c;3;8;74;5;5;oI16;a,b/a,c/a,z2,x3;-;O2Pb1Pd1;O2Pb1Pd1 (ICSD #2277)");
    vproto.push_back("AB3C_oI20_74_e_ef_c;3;10;74;5;6;oI20;a,b/a,c/a,z2,z3,x4;-;Ba1O3Pb1;Ba1O3Pb1 (ICSD #15933)");
    vproto.push_back("A2B2C_oI20_74_ac_i_e;3;10;74;5;6;oI20;a,b/a,c/a,z3,x4,z4;-;Cu2O2Tl1;Cu2O2Tl1 (ICSD #36535)");
    vproto.push_back("ABC4_oI24_74_e_e_hi;3;12;74;5;9;oI24;a,b/a,c/a,z1,z2,y3,z3,x4,z4;-;Cr1Dy1O4;Cr1Dy1O4 (ICSD #93789)");
    vproto.push_back("A2BC4_tI28_79_c_2a_2c;3;14;79;5;13;tI28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Bi2Cu1O4;Bi2Cu1O4 (ICSD #12104)");
    vproto.push_back("ABC_tI24_82_g_g_g;3;12;82;5;11;tI24;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;Ag1K1O1;Ag1K1O1 (ICSD #25744)");
    vproto.push_back("AB5C_tP14_85_a_cg_c;3;14;85;5;7;tP14;a,c/a,z2,z3,x4,y4,z4;-;Mo1O5V1;Mo1O5V1 (ICSD #27315)");
    vproto.push_back("AB4C4_tI18_87_a_h_h;3;9;87;5;6;tI18;a,c/a,x2,y2,x3,y3;-;Ir1Na4O4;Ir1Na4O4 (ICSD #67826)");
    vproto.push_back("ABC_tI24_87_h_h_h;3;12;87;5;8;tI24;a,c/a,x1,y1,x2,y2,x3,y3;-/-;Cu1Na1O1/Cu1O1Rb1;Cu1Na1O1 (ICSD #15099)/Cu1O1Rb1 (ICSD #15100)");
    vproto.push_back("AB8C4_tI26_87_a_2h_h;3;13;87;5;8;tI26;a,c/a,x2,y2,x3,y3,x4,y4;-/-;K1O8Ru4/K1O8V4;K1O8Ru4 (ICSD #1562)/K2O16V8 (ICSD #100596)");
    vproto.push_back("A2B5C_tI32_87_h_di_e;3;16;87;5;8;tI32;a,c/a,z2,x3,y3,x4,y4,z4;-;La4O10Re2;La4O10Re2 (ICSD #81)");
    vproto.push_back("A4BC_tI24_88_f_a_b;3;12;88;5;5;tI24;a,c/a,x3,y3,z3;-;O4Pb1W1;O4Pb1W1 (ICSD #16189)");
    vproto.push_back("AB4C_tI24_88_a_f_b;3;12;88;5;5;tI24;a,c/a,x3,y3,z3;-/-/-;Mo1O4Pb1/Ba1O4W1/Ca1Mo1O4;Mo1O4Pb1 (ICSD #39137)/Ba1O4W1 (ICSD #52384)/Ca1Mo1O4 (ICSD #77334)");
    vproto.push_back("AB4C_tI24_88_b_f_a;3;12;88;5;5;tI24;a,c/a,x3,y3,z3;-;Ba1O4W1;Ba1O4W1 (ICSD #169094)");
    vproto.push_back("A4B2C_tI28_88_f_d_a;3;14;88;5;5;tI28;a,c/a,x3,y3,z3;-;O4Pd2Pr1;O4Pd2Pr1 (ICSD #78843)");
    vproto.push_back("A2BC4_tI28_88_d_a_f;3;14;88;5;5;tI28;a,c/a,x3,y3,z3;-;Au2Ba1O4;Au2Ba1O4 (ICSD #80327)");
    vproto.push_back("ABC2_tP16_92_a_a_b;3;16;92;5;7;tP16;a,c/a,x1,x2,x3,y3,z3;-/-/-;Al1Li1O2/Fe1Na1O2/Al1K1O2;Al1Li1O2 (ICSD #23815)/Fe1Na1O2 (ICSD #33768)/Al1K1O2 (ICSD #151883)");
    vproto.push_back("A3BC_tP5_99_ac_b_a;3;5;99;5;6;tP5;a,c/a,z1,z2,z3,z4;-;O3Pb1V1;O3Pb1V1 (ICSD #152276)");
    vproto.push_back("ABC3_tP5_99_b_a_ac;3;5;99;5;6;tP5;a,c/a,z1,z2,z3,z4;-;Bi1Fe1O3;Bi1Fe1O3 (ICSD #188467)");
    vproto.push_back("A6B5C_tI24_107_bc_ad_a;3;12;107;5;9;tI24;a,c/a,z1,z2,z3,x4,z4,x5,z5;-;Na6O5Pb1;Na6O5Pb1 (ICSD #15102)");
    vproto.push_back("A2B4C_tI28_108_c_d_a;3;14;108;5;8;tI28;a,c/a,z1,x2,z2,x3,y3,z3;-;Bi2O4Pd1;Bi2O4Pd1 (ICSD #9622)");
    vproto.push_back("ABC2_tP8_113_c_a_e;3;8;113;5;5;tP8;a,c/a,z2,x3,z3;-;Ca1Fe1O2;Ca1Fe1O2 (ICSD #246244)");
    vproto.push_back("ABC_tI24_119_g_i_i;3;12;119;5;7;tI24;a,c/a,x1,x2,z2,x3,z3;-/-/-/-/-;Cu1Li1O1/Ag1K1O1/Ag1Na1O1/Ag1Cs1O1/Cu1Li1O1;Cu1Li1O1 (ICSD #282)/Ag1K1O1 (ICSD #24818)/Ag1Na1O1 (ICSD #49752)/Ag1Cs1O1 (ICSD #49754)/Cu1Li1O1 (ICSD #49755)");
    vproto.push_back("A3B4C_tI16_121_bd_i_a;3;8;121;5;4;tI16;a,c/a,x4,z4;-/-;K3O4V1/Ag3O4V1;K3O4V1 (ICSD #4138)/Ag3O4V1 (ICSD #417470)");
    vproto.push_back("A2BC6_tI18_121_e_a_di;3;9;121;5;5;tI18;a,c/a,z3,x4,z4;-;La2Mo1O6;La2Mo1O6 (ICSD #25611)");
    vproto.push_back("AB3C8_tI24_121_a_bd_2i;3;12;121;5;6;tI24;a,c/a,x4,z4,x5,z5;-;Cr1K3O8 (CrK3O8);Cr1K3O8 (ICSD #9356) and K3CrO8 (part 3)"); //DX20210428 - added equivalent part 3 prototype (Kr3CrO8, http://aflow.org/prototype-encyclopedia/AB3C8_tI24_121_a_bd_2i.html)
    vproto.push_back("A3BC8_tI24_121_bd_a_2i;3;12;121;5;6;tI24;a,c/a,x4,z4,x5,z5;-;K3Nb1O8;K3Nb1O8 (ICSD #30405)");
    vproto.push_back("ABC2_tI16_122_b_a_d;3;8;122;5;3;tI16;a,c/a,x3;-;Co1K1O2;Co1K1O2 (ICSD #4199)");
    vproto.push_back("A2BC4_tI28_122_d_a_e;3;14;122;5;6;tI28;a,c/a,x2,x3,y3,z3;-;Cr2Cu1O4;Cr2Cu1O4 (ICSD #16708)");
    vproto.push_back("AB3C_tP5_123_c_ae_b;3;5;123;5;2;tP5;a,c/a;-;Ba1O3Ti1;Ba1O3Ti1 (ICSD #27965)");
    vproto.push_back("AB3C_tP5_123_a_be_d;3;5;123;5;2;tP5;a,c/a;-;Ba1O3Ti1;Ba1O3Ti1 (ICSD #109327)");
    vproto.push_back("ABC3_tP5_123_c_b_ae;3;5;123;5;2;tP5;a,c/a;-;Ca1Mn1O3;Ca1Mn1O3 (ICSD #168904)");
    vproto.push_back("AB4C6_tP11_123_a_eh_bci;3;11;123;5;4;tP11;a,c/a,z5,z6;-;Ba1Nb4O6;Ba1Nb4O6 (ICSD #42006)");
    vproto.push_back("A5B9C2_tP16_123_bfg_cegi_h;3;16;123;5;6;tP16;a,c/a,z5,z6,z7,z8;-;Nb5O9Sr2;Nb5O9Sr2 (ICSD #42007)");
    vproto.push_back("AB3C_tP10_127_c_bg_a;3;10;127;5;3;tP10;a,c/a,x4;-;Na1O3Ta1;Na1O3Ta1 (ICSD #23322)");
    vproto.push_back("AB4C2_tP14_131_c_n_j;3;14;131;5;4;tP14;a,c/a,x2,x3;-;Ca1O4Pt2;Ca1O4Pt2 (ICSD #6159)");
    vproto.push_back("A3BC2_tP12_136_ad_b_f;3;12;136;5;3;tP12;a,c/a,x4;-;K3Ni1O2;K3Ni1O2 (ICSD #262579)");
    vproto.push_back("AB2C2_tI10_139_a_e_d;3;5;139;5;3;tI10;a,c/a,z3;-;Bi1Ce2O2 (Co2S2Tl);Bi1Ce2O2 (ICSD #9099) and Co2S2Tl (part 3)");
    vproto.push_back("AB2C2_tI10_139_a_e_e;3;5;139;5;4;tI10;a,c/a,z2,z3;-;Hg1Na2O2;Hg1Na2O2 (ICSD #25511)");
    vproto.push_back("AB2C4_tI14_139_a_e_cd;3;7;139;5;3;tI14;a,c/a,z4;-/-/-;Cu1Nd2O4/Cu1La2O4/Cu1Nd2O4;Cu1Nd2O4 (ICSD #4203)/Cu1La2O4 (ICSD #261660)/Cu1Nd2O4 (ICSD #86753)");
    vproto.push_back("A2B4C_tI14_139_e_ce_a;3;7;139;5;4;tI14;a,c/a,z3,z4;-;Ba2O4Pb1;Ba2O4Pb1 (ICSD #27113)");
    vproto.push_back("A3B2C6_tI22_139_ae_e_dg;3;11;139;5;5;tI22;a,c/a,z3,z4,z5;-;La3Ni2O6;La3Ni2O6 (ICSD #249209)");
    vproto.push_back("ABC_tI24_139_h_j_i;3;12;139;5;5;tI24;a,c/a,x1,x2,x3;-/-/-/-;Ag1Na1O1/Ag1K1O1/Cu1Li1O1/Ag1Cs1O1;Ag1Na1O1 (ICSD #40153)/Ag1K1O1 (ICSD #40154)/Cu1Li1O1 (ICSD #40156)/Ag1Cs1O1 (ICSD #40160)");
    vproto.push_back("A7B3C2_tI24_139_aeg_be_e;3;12;139;5;6;tI24;a,c/a,z3,z4,z5,z6;-;O7Sr3Ti2 (O7Sr3Ti2);O7Sr3Ti2 (ICSD #20294) and O7Sr3Ti2 (part 3)"); //DX20210428 - added equivalent part 3 prototype (O7Sr3Ti2, http://aflow.org/prototype-encyclopedia/A7B3C2_tI24_139_aeg_be_e.html)
    vproto.push_back("ABC_tI24_139_h_i_j;3;12;139;5;5;tI24;a,c/a,x1,x2,x3;-/-;Cu1O1Rb1/Au1O1Rb1;Cu1O1Rb1 (ICSD #40159)/Au1O1Rb1 (ICSD #409552)");
    vproto.push_back("A4B3C8_tI30_139_2e_ae_cdg;3;15;139;5;6;tI30;a,c/a,z4,z5,z6,z7;-;Nd4Ni3O8;Nd4Ni3O8 (ICSD #51097)");
    vproto.push_back("A3BC_tI20_140_ah_b_c;3;10;140;5;3;tI20;a,c/a,x4;-;O3Sr1Zr1;O3Sr1Zr1 (ICSD #1522)");
    vproto.push_back("AB3C_tI20_140_c_ah_b;3;10;140;5;3;tI20;a,c/a,x4;-;Hf1O3Sr1;Hf1O3Sr1 (ICSD #164621)");
    vproto.push_back("ABC2_tI16_141_b_a_e;3;8;141;5;3;tI16;a,c/a,z3;-;Gd1Na1O2;Gd1Na1O2 (ICSD #22040)");
    vproto.push_back("AB2C2_tI20_141_a_d_e;3;10;141;5;3;tI20;a,c/a,z3;-;Ba1Cu2O2;Ba1Cu2O2 (ICSD #9456)");
    vproto.push_back("AB4C_tI24_141_a_h_b;3;12;141;5;4;tI24;a,c/a,y3,z3;-;Cs1O4Re1;Cs1O4Re1 (ICSD #72817)");
    vproto.push_back("A2B2C3_tI28_141_d_c_ae;3;14;141;5;3;tI28;a,c/a,z4;-;Ag2Cu2O3;Ag2Cu2O3 (ICSD #51502)");
    vproto.push_back("A2B4C_tI28_141_d_h_a;3;14;141;5;4;tI28;a,c/a,y3,z3;-;Mn2O4Zn1;Mn2O4Zn1 (ICSD #15305)");
    vproto.push_back("AB2C4_tI28_141_a_d_h;3;14;141;5;4;tI28;a,c/a,y3,z3;-;Cd1Mn2O4;Cd1Mn2O4 (ICSD #24258)");
    vproto.push_back("A2BC4_tI28_141_d_a_h;3;14;141;5;4;tI28;a,c/a,y3,z3;-/-/-;Cr2Ni1O4/Cr2Fe1O4/Cr2Cu1O4;Cr2Ni1O4 (ICSD #37023)/Cr2Fe1O4 (ICSD #183968)/Cr2Cu1O4 (ICSD #246080)");
    vproto.push_back("ABC2_tI32_141_c_d_h;3;16;141;5;4;tI32;a,c/a,y3,z3;-;Li1Mn1O2;Li1Mn1O2 (ICSD #40486)");
    vproto.push_back("AB2C_tI32_141_d_h_c;3;16;141;5;4;tI32;a,c/a,y3,z3;-;Li1O2Ti1;Li1O2Ti1 (ICSD #164158)");
    vproto.push_back("ABC3_hR5_146_a_a_b;3;5;146;5;7;hR5;a,c/a,x1,x2,x3,y3,z3;-;Bi1Fe1O3;Bi1Fe1O3 (ICSD #162264)");
    vproto.push_back("A3B9C4_hR16_146_3a_3b_4a;3;16;146;5;18;hR16;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-/-;Ba3O9Yb4/Ba3O9Y4;Ba3O9Yb4 (ICSD #33239)/Ba3O9Y4 (ICSD #87118)");
    vproto.push_back("A3B4C9_hR16_146_3a_4a_3b;3;16;146;5;18;hR16;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Ba3Ho4O9;Ba3Ho4O9 (ICSD #33807)");
    vproto.push_back("AB8C2_hP11_147_a_dg_d;3;11;147;5;7;hP11;a,c/a,z2,z3,x4,y4,z4;-/-;Mn1O8Re2/Ni1O8Re2;Mn1O8Re2 (ICSD #51014)/Ni1O8Re2 (ICSD #51016)");
    vproto.push_back("A8B6C_hR15_148_cf_f_a;3;15;148;5;9;hR15;a,c/a,x2,x3,y3,z3,x4,y4,z4;-;Li8O6Sn1 (Li7O6Ta);Li8O6Sn1 (ICSD #1180) and Li7O6Ta (partially occupied Li, part 3)"); //DX20210428 - added equivalent part 3 prototype (Li7O6Ta, partially occupied Li, http://aflow.org/prototype-encyclopedia/A8B6C_hR15_148_cf_f_a.html)
    vproto.push_back("AB2C_hP12_152_b_c_a;3;12;152;5;7;hP12;a,c/a,x1,x2,x3,y3,z3;-;Ba1O2Zn1;Ba1O2Zn1 (ICSD #25559)");
    vproto.push_back("AB2C_hP12_152_a_c_b;3;12;152;5;7;hP12;a,c/a,x1,x2,x3,y3,z3;-;Ba1O2Zn1;Ba1O2Zn1 (ICSD #25812)");
    vproto.push_back("AB2C_hP12_154_a_c_b;3;12;154;5;7;hP12;a,c/a,x1,x2,x3,y3,z3;-;Hg1O2Sr1;Hg1O2Sr1 (ICSD #69739)");
    vproto.push_back("AB3C_hR5_160_a_b_a;3;5;160;5;6;hR5;a,c/a,x1,x2,x3,z3;-;Ba1O3Ti1;Ba1O3Ti1 (ICSD #6102)");
    vproto.push_back("A2B3C2_hR7_160_2a_b_2a;3;7;160;5;8;hR7;a,c/a,x1,x2,x3,x4,x5,z5;-;K2O3Sn2;K2O3Sn2 (ICSD #2216)");
    vproto.push_back("AB3C_hR10_161_a_b_a;3;10;161;5;7;hR10;a,c/a,x1,x2,x3,y3,z3;-;Ni1O3Pb1;Ni1O3Pb1 (ICSD #187684)");
    vproto.push_back("A2B3C_hR12_161_2a_b_a;3;12;161;5;8;hR12;a,c/a,x1,x2,x3,x4,y4,z4;-;Li2O3Re1;Li2O3Re1 (ICSD #200999)");
    vproto.push_back("A5B6C2_hP13_162_ef_k_d;3;13;162;5;5;hP13;a,c/a,z2,x4,z4;-;Ag5O6Pb2;Ag5O6Pb2 (ICSD #40058)");
    vproto.push_back("A4B7C2_hP13_162_de_ak_c;3;13;162;5;5;hP13;a,c/a,z4,x5,z5;-;In4O7Rb2;In4O7Rb2 (ICSD #6321)");
    vproto.push_back("AB3C2_hP6_164_b_e_d;3;6;164;5;3;hP6;a,c/a,z2;-;Cs1Cu3O2;Cs1Cu3O2 (ICSD #413342)");
    vproto.push_back("A2B8C_hP11_164_d_di_a;3;11;164;5;6;hP11;a,c/a,z2,z3,x4,z4;-;Mo2O8Zr1;Mo2O8Zr1 (ICSD #59999)");
    vproto.push_back("A7B4C2_hP13_164_ai_2d_c;3;13;164;5;7;hP13;a,c/a,z2,z3,z4,x5,z5;-;O7Tl4V2;O7Tl4V2 (ICSD #72809)");
    vproto.push_back("A3B8C3_hP14_164_ad_di_bd;3;14;164;5;7;hP14;a,c/a,z3,z4,z5,x6,z6;-;K3O8V3;K3O8V3 (ICSD #100782)");
    vproto.push_back("A4B2C_hP14_164_di_abd_d;3;14;164;5;7;hP14;a,c/a,z3,z4,z5,x6,z6;-;O4Tl2W1;O4Tl2W1 (ICSD #8212)");
    vproto.push_back("ABC2_hR4_166_b_a_c;3;4;166;5;3;hR4;a,c/a,x3;-/-/-/-/-;Al1Cu1O2/Mo1Na1O2/C1Ca1N2/C1Ca1N2/C1Mg1N2;Al1Cu1O2 (ICSD #25593)/Mo1Na1O2 (ICSD #166514)/C1Ca1N2 (ICSD #25763)/C1Ca1N2 (ICSD #31100)/C1Mg1N2 (ICSD #75039)"); //DX20210221 - added three carbo-nitrides
    vproto.push_back("A2BC_hR4_166_c_a_b;3;4;166;5;3;hR4;a,c/a,x3;-;O2Rb1Sc1;O2Rb1Sc1 (ICSD #31959)");
    vproto.push_back("ABC3_hR5_166_b_a_d;3;5;166;5;2;hR5;a,c/a;-;Bi1Fe1O3;Bi1Fe1O3 (ICSD #20288)");
    vproto.push_back("A2BC2_hR5_166_c_a_c;3;5;166;5;4;hR5;a,c/a,x2,x3;-;Ag2Ni1O2;Ag2Ni1O2 (ICSD #160574)");
    vproto.push_back("A2B4C_hR7_166_c_2c_a;3;7;166;5;5;hR7;a,c/a,x2,x3,x4;-/-;Fe2O4Yb1/Fe2O4Y1;Fe2O4Yb1 (ICSD #4192)/Fe2O4Y1 (ICSD #67701)");
    vproto.push_back("A2B3C2_hR7_166_ab_d_c;3;7;166;5;3;hR7;a,c/a,x3;-;K2O3Sn2 (Ni3Pb2S2);K2O3Sn2 (ICSD #15511) and Ni3Pb2S2 (shandite, part 3)"); //DX20210428 - added equivalent part 3 prototype (Ni3Pb2S2, shandite, http://aflow.org/prototype-encyclopedia/A3B2C2_hR7_166_d_ab_c.html)
    vproto.push_back("A2BC4_hR7_166_c_a_2c;3;7;166;5;5;hR7;a,c/a,x2,x3,x4;-/-;Fe2In1O4/Bi2MnTe4;Fe2In1O4 (ICSD #157323)/MnBi2Te4"); //DX20210427 - added MnBi2Te4
    vproto.push_back("AB2C_hR8_166_c_ad_c;3;8;166;5;4;hR8;a,c/a,x2,x3;-;Al1O2Tl1;Al1O2Tl1 (ICSD #29010)");
    vproto.push_back("AB3C2_hR12_166_c_h_ae;3;12;166;5;5;hR12;a,c/a,x2,x4,z4;-;K1O3Pd2;K1O3Pd2 (ICSD #248051)");
    vproto.push_back("AB4C8_hR13_166_a_bd_ch;3;13;166;5;5;hR13;a,c/a,x3,x5,z5;-/-;Ba1Ni4O8/Ca1Ni4O8;Ba1Ni4O8 (ICSD #20898)/Ca1Ni4O8 (ICSD #40470)");
    vproto.push_back("A3B2C8_hR13_166_ac_c_ch;3;13;166;5;7;hR13;a,c/a,x2,x3,x4,x5,z5;-;Ba3Cr2O8;Ba3Cr2O8 (ICSD #9457)");
    vproto.push_back("A2B8C3_hR13_166_c_ch_e;3;13;166;5;6;hR13;a,c/a,x1,x2,x4,z4;-;Ca2O8Pt3;Ca2O8Pt3 (ICSD #65412)");
    vproto.push_back("A4B8C_hR13_166_bd_ch_a;3;13;166;5;5;hR13;a,c/a,x3,x5,z5;-;Ni4O8Sr1;Ni4O8Sr1 (ICSD #40469)");
    vproto.push_back("A8B3C2_hR13_166_ch_ac_c;3;13;166;5;7;hR13;a,c/a,x2,x3,x4,x5,z5;-;O8Pb3V2;O8Pb3V2 (ICSD #27651)");
    vproto.push_back("AB4C2_hR14_166_c_ch_ad;3;14;166;5;6;hR14;a,c/a,x2,x3,x5,z5;-;Al1O4V2;Al1O4V2 (ICSD #151457)");
    vproto.push_back("AB3C_hR15_166_bc_dh_ac;3;15;166;5;6;hR15;a,c/a,x3,x4,x6,z6;-;Ba1O3Ru1;Ba1O3Ru1 (ICSD #10253)");
    vproto.push_back("A3BC4_hR16_167_e_a_be;3;16;167;5;4;hR16;a,c/a,x3,x4;-;Ba3Ni1O4;Ba3Ni1O4 (ICSD #30662)");
    vproto.push_back("AB2C4_hP14_182_b_f_cg;3;14;182;5;4;hP14;a,c/a,z3,x4;-;Ba1Fe2O4;Ba1Fe2O4 (ICSD #2769)");
    vproto.push_back("ABC3_hP10_186_b_a_c;3;10;186;5;6;hP10;a,c/a,z1,z2,x3,z3;-;Ba1Ni1O3;Ba1Ni1O3 (ICSD #15761)");
    vproto.push_back("AB3C_hP10_186_a_a2b_b;3;10;186;5;7;hP10;a,c/a,z1,z2,z3,z4,z5;-;La1O3Tl1;La1O3Tl1 (ICSD #200088)");
    vproto.push_back("AB2C4_hP14_186_b_ab_bc;3;14;186;5;8;hP14;a,c/a,z1,z2,z3,z4,x5,z5;-;Mn1Na2O4;Mn1Na2O4 (ICSD #39504)");
    vproto.push_back("AB2C_hP8_187_bc_gh_i;3;8;187;5;5;hP8;a,c/a,z3,z4,z5;-;Cs1O2Y1;Cs1O2Y1 (ICSD #49652)");
    vproto.push_back("ABC3_hP10_187_ad_i_jk;3;10;187;5;5;hP10;a,c/a,z3,x4,x5;-;Ba1Co1O3;Ba1Co1O3 (ICSD #88670)");
    vproto.push_back("A3B4C_hP16_190_af_cg_d;3;16;190;5;4;hP16;a,c/a,z4,x5;-;Be3O4Sr1;Be3O4Sr1 (ICSD #26179)");
    vproto.push_back("A2B3C_hP6_191_d_f_a;3;6;191;5;2;hP6;a,c/a;-;Mn2O3Ta1;Mn2O3Ta1 (ICSD #15995)");
    vproto.push_back("A2BC_hP8_194_f_c_a;3;8;194;5;3;hP8;a,c/a,z3;-/-;O2Rb1Sc1/C2Cr1Sc1;O2Rb1Sc1 (ICSD #1270)/C2Cr1Sc1 (ICSD #80373)"); //DX20210121 - added metal-carbide
    vproto.push_back("AB2C_hP8_194_c_f_a;3;8;194;5;3;hP8;a,c/a,z3;-/-/-;Cu1O2Y1/Cu1O2Sc1/C1N2Ni1;Cu1O2Y1 (ICSD #35580)/Cu1O2Sc1 (ICSD #60847)/C1N2Ni1 (ICSD #249388)"); //DX20210221 - added carbo-nitride
    vproto.push_back("AB3C_hP10_194_c_bf_a;3;10;194;5;3;hP10;a,c/a,z4;-/-;Al1O3Y1/Mn1O3Y1;Al1O3Y1 (ICSD #27100)/Mn1O3Y1 (ICSD #73361)");
    vproto.push_back("ABC3_hP10_194_c_a_bf;3;10;194;5;3;hP10;a,c/a,z4;-;Ga1In1O3;Ga1In1O3 (ICSD #30339)");
    vproto.push_back("ABC3_hP10_194_a_c_bf;3;10;194;5;3;hP10;a,c/a,z4;-;In1Mn1O3;In1Mn1O3 (ICSD #67671)");
    vproto.push_back("ABC2_hP16_194_ac_f_ef;3;16;194;5;5;hP16;a,c/a,z3,z4,z5;-;Ag1Co1O2;Ag1Co1O2 (ICSD #261608)");
    vproto.push_back("A2B3C2_cI28_199_a_b_a;3;14;199;5;4;cI28;a,x1,x2,x3;-;K2O3Pb2;K2O3Pb2 (ICSD #1412)");
    vproto.push_back("A4B2C7_cF52_216_e_bc_ag;3;13;216;5;3;cF52;a,x4,x5;-;Al4Cu2O7;Al4Cu2O7 (ICSD #100355)");
    vproto.push_back("A4B3C_cP16_223_e_c_a;3;16;223;5;1;cP16;a;-;O4Pd3Sr1;O4Pd3Sr1 (ICSD #16537)");
    vproto.push_back("ABC3_cF40_225_c_ab_e;3;10;225;5;2;cF40;a,x4;-;Ba2Bi2O6;Ba2Bi2O6 (ICSD #68714)");
    vproto.push_back("A3B6C_cF40_225_bc_e_a;3;10;225;5;2;cF40;a,x4;-/-;Eu3O6Ta1/Ba3O6W1;Eu3O6Ta1 (ICSD #4148)/Ba3O6W1 (ICSD #76437)");
    vproto.push_back("A6BC8_cF60_225_d_a_ce;3;15;225;5;2;cF60;a,x4;-;Mg6Mn1O8;Mg6Mn1O8 (ICSD #24710)");
    vproto.push_back("A6B8C_cF60_225_d_f_a;3;15;225;5;2;cF60;a,x3;-;Cu6O8Pb1;Cu6O8Pb1 (ICSD #280596)");
    vproto.push_back("A4B3C_cF64_225_f_d_ab;3;16;225;5;2;cF64;a,x4;-;O4Pd3Tl1;O4Pd3Tl1 (ICSD #2275)");
    vproto.push_back("AB2C_cF64_227_c_e_d;3;16;227;5;2;cF64;a,x3;-;Li1O2Ti1;Li1O2Ti1 (ICSD #48128)");
    // -------------------------------------------------------------------------
    // metal-carbide prototypes (from DX) //DX2020120
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A2B_aP12_2_4i_2i;2;12;2;4;24;aP12;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;C2Ca1;C2Ca1 (ICSD #66663)");
    vproto.push_back("AB4_mP10_11_e_4e;2;10;11;4;14;mP10;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;-;C1Fe4;C1Fe4 (ICSD #187143)");
    vproto.push_back("A5B6_mC22_12_agh_ij;2;11;12;4;11;mC22;a,b/a,c/a,beta,y2,y3,x4,z4,x5,y5,z5;-;C5Nb6;C5Nb6 (ICSD #20695)");
    vproto.push_back("A4B_mP10_14_2e_a;2;10;14;4;10;mP10;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3;-;C4Ir1;C4Ir1 (ICSD #181493)");
    vproto.push_back("A2B_mC12_15_f_e;2;6;15;4;8;mC12;a,b/a,c/a,beta,y1,x2,y2,z2;-;C2Ca1;C2Ca1 (ICSD #54184)");
    vproto.push_back("A3B5_oP16_55_ah_d2g;2;16;55;4;9;oP16;a,b/a,c/a,x3,y3,x4,y4,x5,y5;-;C3Ir5;C3Ir5 (ICSD #181485)");
    vproto.push_back("A3B2_oP10_58_ag_g;2;10;58;4;7;oP10;a,b/a,c/a,x2,y2,x3,y3;-;C3Mg2;C3Mg2 (ICSD #71941)");
    vproto.push_back("AB_oP16_62_2c_2c;2;16;62;4;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;-;C2Rb2;C2Rb2 (ICSD #51529)");
    vproto.push_back("A3B_oP16_62_3c_c;2;16;62;4;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;-/D0_{16};C3Ir1/I3(NH4);C3Ir1 (ICSD #181492)/NH4I3 (D016)"); //DX20210427 - added NH4I3 (D016) from part 3
    vproto.push_back("A4B_oC20_68_i_a;2;10;68;4;6;oC20;a,b/a,c/a,x2,y2,z2;-;C4Ir1;C4Ir1 (ICSD #181496)");
    vproto.push_back("A2B_tP3_123_h_a;2;3;123;4;3;tP3;a,c/a,z2;-;C2Re1;C2Re1 (ICSD #184664)");
    vproto.push_back("AB_tI32_142_f_e;2;16;142;4;4;tI32;a,c/a,x1,x2;-/-;C2Na2/C2K2;C2Na2 (ICSD #28066)/C2K2 (ICSD #36142)");
    vproto.push_back("A4B_hP15_152_2ac_b;2;15;152;4;8;hP15;a,c/a,x1,x2,x3,x4,y4,z4;-;C4Ir1;C4Ir1 (ICSD #181495)");
    vproto.push_back("A4B3_hR7_160_4a_3a;2;7;160;4;9;hR7;a,c/a,x1,x2,x3,x4,x5,x6,x7;-;Al4C3;Al4C3 (ICSD #14397)");
    vproto.push_back("A5B8_hR13_160_2ab_2a2b;2;13;160;4;12;hR13;a,c/a,x1,x2,x3,x4,x5,z5,x6,z6,x7,z7;-;C5Ti8;C5Ti8 (ICSD #618924)");
    vproto.push_back("A5B8_hR13_166_abd_ch;2;13;166;4;5;hR13;a,c/a,x3,x5,z5;-;C5Ti8;C5Ti8 (ICSD #20822)");
    vproto.push_back("AB2_hP3_187_a_i;2;3;187;4;3;hP3;a,c/a,z2;-;C1Mo2;C1Mo2 (ICSD #43669)");
    vproto.push_back("AB2_hP3_191_a_c;2;6;191;4;2;hP3;a,c/a;-;C1Fe2;C1Fe2 (ICSD #162103)");
    vproto.push_back("A8B_hP9_191_cl_b;2;9;191;4;3;hP9;a,c/a,x3;-;C8Cs1;C8Cs1 (ICSD #74641)");
    vproto.push_back("A6B_hP14_194_i_c;2;14;194;4;3;hP14;a,c/a,x2;-/-;C6Eu1/C6Yb1;C6Eu1 (ICSD #169041)/C6Yb1 (ICSD #601565)");
    vproto.push_back("A3B4_cI28_220_a_c;2;14;220;4;2;cI28;a,x2;-;C3Sc4;C3Sc4 (ICSD #42760)");
    vproto.push_back("A2B_cP12_224_e_b;2;12;224;4;2;cP12;a,x2;-;C2Ca1;C2Ca1 (ICSD #31092)");
    vproto.push_back("AB_cF64_227_cd_e;2;16;227;4;2;cF64;a,x3;-;C1Ti1;C1Ti1 (ICSD #618950)");
    // -------------------------------------------------------------------------
    // ternaries
    vproto.push_back("A2BC_mC16_9_2a_a_a;3;8;9;5;16;mC16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-/-;C2Co1Nd1/C2Co1Nd1/C2Co1Pr1;C2Co1Nd1 (ICSD #55570)/C2Co1Nd1 (ICSD #67375)/C2Co1Pr1 (ICSD #658525)");
    vproto.push_back("A5B4C3_mC24_12_b2i_2i_ai;3;12;12;5;14;mC24;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;C5Ca4Ni3;C5Ca4Ni3 (ICSD #71440)");
    vproto.push_back("A4BC3_mC32_12_2j_i_ghi;3;16;12;5;16;mC32;a,b/a,c/a,beta,y1,y2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6;-;C4Co1Sc3;C4Co1Sc3 (ICSD #236391)");
    vproto.push_back("A4BC3_mC32_12_2j_i_bc2i;3;16;12;5;16;mC32;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;C4Ru1Sc3;C4Ru1Sc3 (ICSD #420074)");
    vproto.push_back("A5B2C4_oP11_25_c2g_g_adh;3;11;25;5;14;oP11;a,b/a,c/a,z1,z2,z3,y4,z4,y5,z5,y6,z6,y7,z7;-;C5Ni2Yb4;C5Ni2Yb4 (ICSD #73156)");
    vproto.push_back("A2BC_oC8_38_d_b_a;3;4;38;5;7;oC8;a,b/a,c/a,z1,z2,y3,z3;-;C2Nd1Rh1;C2Nd1Rh1 (ICSD #63064)");
    vproto.push_back("A2BC_oP16_57_2d_d_c;3;16;57;5;10;oP16;a,b/a,c/a,x1,x2,y2,x3,y3,x4,y4;-;C2Cs1Na1;C2Cs1Na1 (ICSD #189824)");
    vproto.push_back("A2BC_oP16_59_2e_ab_e;3;16;59;5;11;oP16;a,b/a,c/a,z1,z2,y3,z3,y4,z4,y5,z5;-;C2Cr1Sc1;C2Cr1Sc1 (ICSD #80372)");
    vproto.push_back("A2BC_oC16_63_g_c_a;3;8;63;5;6;oC16;a,b/a,c/a,y2,x3,y3;-;C2Gd1Ru1;C2Gd1Ru1 (ICSD #80312)");
    vproto.push_back("A2B2C_oC20_63_ac_f_c;3;10;63;5;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;C2Cr2V1;C2Cr2V1 (ICSD #20297)");
    vproto.push_back("A4BC3_oC16_65_p_a_bi;3;8;65;5;6;oC16;a,b/a,c/a,y3,x4,y4;-;C4Os1Sc3;C4Os1Sc3 (ICSD #420076)");
    vproto.push_back("A7B2C5_oC28_65_aegj_i_bq;3;14;65;5;8;oC28;a,b/a,c/a,x4,y5,y6,x7,y7;-;C7Re2Sc5;C7Re2Sc5 (ICSD #72382)");
    vproto.push_back("A4BC3_oI16_71_n_b_af;3;8;71;5;6;oI16;a,b/a,c/a,x3,x4,y4;-;C4Fe1Sc3;C4Fe1Sc3 (ICSD #72863)");
    vproto.push_back("A4B2C_oI28_72_2j_j_a;3;14;72;5;9;oI28;a,b/a,c/a,x2,y2,x3,y3,x4,y4;-;C4Er2Fe1;C4Er2Fe1 (ICSD #42968)");
    vproto.push_back("A2BC_tP4_123_g_a_d;3;4;123;5;3;tP4;a,c/a,z3;-;C2Cu1Rb1;C2Cu1Rb1 (ICSD #391118)");
    vproto.push_back("AB2C_tP4_123_a_g_d;3;4;123;5;3;tP4;a,c/a,z3;-/-;Ag1C2K1/Au1C2Na1;Ag1C2K1 (ICSD #410874)/Au1C2Na1 (ICSD #411254)");
    vproto.push_back("AB2C2_tP5_123_b_e_ac;3;5;123;5;2;tP5;a,c/a;-;C1Co2Mn2;C1Co2Mn2 (ICSD #44353)");
    vproto.push_back("A2BC_tP8_129_2c_a_c;3;8;129;5;5;tP8;a,c/a,z2,z3,z4;-;C2Co1Sc1;C2Co1Sc1 (ICSD #62598)");
    vproto.push_back("AB2C_tP8_131_a_j_f;3;8;131;5;3;tP8;a,c/a,x3;-;Ag1C2Cs1;Ag1C2Cs1 (ICSD #410873)");
    vproto.push_back("A2BC_tP8_131_j_a_f;3;8;131;5;3;tP8;a,c/a,x3;-;C2Cu1K1;C2Cu1K1 (ICSD #412038)");
    vproto.push_back("A4B5C2_hR11_160_4a_5a_2a;3;11;160;5;13;hR11;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;-;Al4C5Zr2;Al4C5Zr2 (ICSD #173676)");
    vproto.push_back("A4B6C3_hR13_160_4a_6a_3a;3;13;160;5;15;hR13;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13;-;Al4C6Zr3;Al4C6Zr3 (ICSD #173677)");
    vproto.push_back("A2B2C_hP5_164_c_d_a;3;5;164;5;4;hP5;a,c/a,z2,z3;-;C2Na2Pd1;C2Na2Pd1 (ICSD #50172)");
    vproto.push_back("A2BC2_hP5_164_c_a_d;3;5;164;5;4;hP5;a,c/a,z2,z3;-;C2Pd1Rb2;C2Pd1Rb2 (ICSD #94394)");
    vproto.push_back("A4B5C2_hR11_166_2c_a2c_c;3;11;166;5;7;hR11;a,c/a,x2,x3,x4,x5,x6;-;Al4C5Hf2;Al4C5Hf2 (ICSD #161586)");
    vproto.push_back("A4B6C3_hR13_166_2c_3c_ac;3;13;166;5;8;hR13;a,c/a,x2,x3,x4,x5,x6,x7;-;Al4C6Hf3;Al4C6Hf3 (ICSD #161585)");
    vproto.push_back("A3B3C_hP14_186_3b_3b_a;3;14;186;5;9;hP14;a,c/a,z1,z2,z3,z4,z5,z6,z7;-;Al3C3Sc1;Al3C3Sc1 (ICSD #43477)");
    vproto.push_back("AB2C_hP4_187_a_g_d;3;4;187;5;3;hP4;a,c/a,z3;-;Ag1C2Li1;Ag1C2Li1 (ICSD #410868)");
    vproto.push_back("AB2C3_hP12_194_b_f_af;3;12;194;5;4;hP12;a,c/a,z3,z4;-;Al1C2Ti3;Al1C2Ti3 (ICSD #93503)");
    vproto.push_back("A2BC3_hP12_194_f_b_af;3;12;194;5;4;hP12;a,c/a,z3,z4;-;C2Sn1Ti3;C2Sn1Ti3 (ICSD #160572)");
    vproto.push_back("A3B3C_hP14_194_cf_df_a;3;14;194;5;4;hP14;a,c/a,z4,z5;-;Al3C3Sc1;Al3C3Sc1 (ICSD #62308)");
    // -------------------------------------------------------------------------
    // carbo-nitride prototypes (from DX) //DX2020219
    // -------------------------------------------------------------------------
    // binaries
    vproto.push_back("A11B4_oP15_16_acgqtu_u;2;15;16;4;11;oP15;a,b/a,c/a,z4,z5,x6,y6,z6,x7,y7,z7;-;C11N4;C11N4 (ICSD #184897)");
    vproto.push_back("AB2_oC12_36_a_2a;2;6;36;4;9;oC12;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-/-;C1N2/MoP2;C1N2 (ICSD #247680)/MoP2"); //DX20210427 - added MoP2 from part 3
    vproto.push_back("A11B4_tP15_111_abcmn_n;2;15;111;4;7;tP15;a,c/a,z4,x5,z5,x6,z6;-;C11N4;C11N4 (ICSD #184896)");
    vproto.push_back("AB2_tP6_113_a_e;2;6;113;4;4;tP6;a,c/a,x2,z2;-;C1N2;C1N2 (ICSD #247676)");
    vproto.push_back("AB2_tI6_119_a_f;2;3;119;4;3;tI6;a,c/a,z2;-;C1N2;C1N2 (ICSD #247678)");
    vproto.push_back("AB2_tI24_122_d_e;2;12;122;4;6;tI24;a,c/a,x1,x2,y2,z2;-;C1N2;C1N2 (ICSD #247677)");
    vproto.push_back("A3B4_hR7_160_b_ab;2;7;160;4;7;hR7;a,c/a,x1,x2,z2,x3,z3;-;C3N4;C3N4 (ICSD #41952)");
    vproto.push_back("AB2_hP6_164_c_2d;2;6;164;4;5;hP6;a,c/a,z1,z2,z3;-;C1N2;C1N2 (ICSD #247679)");
    vproto.push_back("A3B4_hP14_176_h_ch;2;14;176;4;6;hP14;a,c/a,x2,y2,x3,y3;-;C3N4 (N4Si3);C3N4 (ICSD #41950) and beta-Si3N4 (part 3)"); //DX20210428 - added equivalent part 3 prototype (beta-Si3N4, http://aflow.org/prototype-encyclopedia/A4B3_hP14_176_ch_h.html)
    vproto.push_back("A3B4_hP14_187_jk_adjk;2;14;187;4;6;hP14;a,c/a,x3,x4,x5,x6;-;C3N4;C3N4 (ICSD #83265)");
    vproto.push_back("A3B4_cP7_215_c_e;2;7;215;4;2;cP7;a,x2;-;C3N4;C3N4 (ICSD #41951)");
    vproto.push_back("A3B4_cI28_220_b_c;2;14;220;4;2;cI28;a,x2;-;C3N4;C3N4 (ICSD #83263)");
    // -------------------------------------------------------------------------
    // ternaries
    vproto.push_back("AB2C2_aP15_2_ai_3i_3i;3;15;2;5;27;aP15;a,b/a,c/a,alpha,beta,gamma,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;C1N2Tl2;C1N2Tl2 (ICSD #417297)");
    vproto.push_back("A4B3C_aP16_2_4i_3i_i;3;16;2;5;30;aP16;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;C4N3Na1;C4N3Na1 (ICSD #31926)");
    vproto.push_back("ABC_mC6_8_a_a_a;3;3;8;5;10;mC6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;C1K1N1;C1K1N1 (ICSD #27350)");
    vproto.push_back("ABC_mC12_9_a_a_a;3;6;9;5;13;mC12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;C1K1N1;C1K1N1 (ICSD #173942)");
    vproto.push_back("AB2C_oP16_33_a_2a_a;3;16;33;5;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;C1N2Pb1;C1N2Pb1 (ICSD #16600)");
    vproto.push_back("AB4C3_oI32_46_b_2bc_bc;3;16;46;5;17;oI32;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,x5,y5,z5,x6,y6,z6;-;Ag1C4N3;Ag1C4N3 (ICSD #43823)");
    vproto.push_back("ABC_oP6_59_a_a_b;3;6;59;5;6;oP6;a,b/a,c/a,z1,z2,z3;-;C1N1Na1;C1N1Na1 (ICSD #77172)");
    vproto.push_back("ABC2_oC16_63_c_a_f;3;8;63;5;6;oC16;a,b/a,c/a,y2,y3,z3;-;C1Cu1N2;C1Cu1N2 (ICSD #161460)");
    vproto.push_back("AB2C_tI32_122_d_e_d;3;16;122;5;7;tI32;a,c/a,x1,x2,x3,y3,z3;-;C1N2Zn1;C1N2Zn1 (ICSD #280523)");
    vproto.push_back("A3B2C6_hR11_155_d_c_f;3;11;155;5;7;hR11;a,c/a,x1,y2,x3,y3,z3;-;C3Lu2N6;C3Lu2N6 (ICSD #240311)");
    vproto.push_back("ABC2_hR4_160_a_a_2a;3;4;160;5;6;hR4;a,c/a,x1,x2,x3,x4;-;C1Cd1N2;C1Cd1N2 (ICSD #95265)");
    vproto.push_back("A6B3C2_hR11_160_6a_3a_2a;3;11;160;5;13;hR11;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;-;Al6C3N2;Al6C3N2 (ICSD #14399)");
    vproto.push_back("A8B3C4_hR15_160_8a_3a_4a;3;15;160;5;17;hR15;a,c/a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15;-;Al8C3N4;Al8C3N4 (ICSD #14401)");
    vproto.push_back("A6B3C2_hR11_166_3c_ac_c;3;11;166;5;7;hR11;a,c/a,x2,x3,x4,x5,x6;-;Al6C3N2;Al6C3N2 (ICSD #41260)");
    vproto.push_back("A8B3C4_hR15_166_4c_ac_2c;3;15;166;5;9;hR15;a,c/a,x2,x3,x4,x5,x6,x7,x8;-;Al8C3N4;Al8C3N4 (ICSD #41261)");
    vproto.push_back("A2BC2_cP10_215_e_ab_e;3;10;215;5;3;cP10;a,x3,x4;-;C2Cd1N2;C2Cd1N2 (ICSD #20748)");
    // -------------------------------------------------------------------------
    // Part 3 //DX20210427
    // -------------------------------------------------------------------------
    vproto.push_back("A5B11CD8E_aP26_1_5a_11a_a_8a_a;5;26;1;7;84;aP26;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26;-;C5H11NaO8S;NaC5H11O8S");
    vproto.push_back("A2B2C5_aP18_2_2i_2i_5i;3;18;2;5;33;aP18;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;B2Co2O5;Co2B2O5");
    vproto.push_back("ABC8D3_aP26_2_i_i_8i_3i;4;26;2;6;45;aP26;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;S6_{8};AlNaO8Si3;Albite (NaAlSi3O8, S68)");
    vproto.push_back("AB3C3_aP28_2_2i_6i_6i;3;28;2;5;48;aP28;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;G5_{1};BH3O3;Boric Acid (H3BO3, G51)");
    vproto.push_back("A3B_aP32_2_12i_4i;2;32;2;4;54;aP32;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;O3W;delta-WO3");
    vproto.push_back("A3B2C10D3_aP36_2_ah2i_2i_10i_3i;4;36;2;6;57;aP36;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19;-;Co3H2O10Se3;Co3(SeO3)3.H2O");
    vproto.push_back("AB10C9D_aP42_2_ae_10i_9i_i;4;42;2;6;66;aP42;a,b/a,c/a,alpha,beta,gamma,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;H4_{10};CuH10O9S;Chalcanthite (CuSO4.5H2O, H410)");
    vproto.push_back("A2B7C2_aP44_2_4i_14i_4i;3;44;2;5;72;aP44;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;-;Ho2O7Si2;alpha-Ho2Si2O7");
    vproto.push_back("A12B2CD12_aP54_2_12i_2i_i_12i;4;54;2;6;87;aP54;a,b/a,c/a,alpha,beta,gamma,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27;-;H12N2NiO12;Ni(NO3)2(H2O)6");
    vproto.push_back("A2B2C5D_mP20_4_2a_2a_5a_a;4;20;4;6;34;mP20;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;H4_{8};H2Li2O5S;Li2SO4.H2O (H48)");
    vproto.push_back("A3B6C_mP20_4_3a_6a_a;3;20;4;5;34;mP20;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Ca3O6U;Ca3UO6");
    vproto.push_back("A11B2C2_mP60_4_22a_4a_4a;3;60;4;5;94;mP60;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30;-;O11P2W2;W2O3(PO4)2");
    vproto.push_back("ABC3_mC10_5_b_a_ac;3;5;5;5;10;mC10;a,b/a,c/a,beta,y1,y2,y3,x4,y4,z4;-;(Ba,Ca)CO3;C2 (Ba,Ca)CO3");
    vproto.push_back("A2B_mC12_5_2c_c;2;6;5;4;13;mC12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;As2Nb;NbAs2");
    vproto.push_back("AB3_mC16_5_c_3c;2;8;5;4;16;mC16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;D0_{15};AlCl3;D015 (AlCl3) (obsolete)");
    vproto.push_back("AB6C18D4E2_mC62_5_a_2b2c_9c_2c_c;5;31;5;7;49;mC62;a,b/a,c/a,beta,y1,y2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;-;CaCu6O18P4Rb2;Rb2CaCu6(PO4)4O2");
    vproto.push_back("A2B2C9D2_mC90_5_ab2c_3c_b13c_3c;4;45;5;6;70;mC90;a,b/a,c/a,beta,y1,y2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24;H4_{7};Ca2H2O9S2;Bassanite [CaSO4(H2O)0.5, H47]");
    vproto.push_back("AB2_mP12_7_2a_4a;2;12;7;4;22;mP12;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;AuTe2;Calaverite (AuTe2)");
    vproto.push_back("A6B2C15D4_mP54_7_6a_2a_15a_4a;4;54;7;6;85;mP54;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27;-;Ca6Na2O15Si4;Na2Ca6Si4O15");
    vproto.push_back("A8B23_mP124_7_16a_46a;2;124;7;4;190;mP124;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36,x37,y37,z37,x38,y38,z38,x39,y39,z39,x40,y40,z40,x41,y41,z41,x42,y42,z42,x43,y43,z43,x44,y44,z44,x45,y45,z45,x46,y46,z46,x47,y47,z47,x48,y48,z48,x49,y49,z49,x50,y50,z50,x51,y51,z51,x52,y52,z52,x53,y53,z53,x54,y54,z54,x55,y55,z55,x56,y56,z56,x57,y57,z57,x58,y58,z58,x59,y59,z59,x60,y60,z60,x61,y61,z61,x62,y62,z62;-;Mo8O23;Low-Temperature Mo8O23");
    vproto.push_back("ABC2_mC8_8_a_a_b;3;4;8;5;11;mC8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,y3,z3;F5_{11};KNO2;F511 (KNO2) (obsolete)");
    vproto.push_back("A13B4_mC102_8_17a11b_8a2b;2;51;8;4;93;mC102;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,z15,x16,z16,x17,z17,x18,z18,x19,z19,x20,z20,x21,z21,x22,z22,x23,z23,x24,z24,x25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36,x37,y37,z37,x38,y38,z38;-;Al13Co4;Monoclinic Co4Al13");
    vproto.push_back("A3B5C4D2_mC56_9_3a_5a_4a_2a;4;28;9;6;46;mC56;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Mg3O5(OH)4Si2;Chrysotile (Mg3Si2O5(OH)4)");
    vproto.push_back("A2B4C9D2_mC68_9_2a_4a_9a_2a;4;34;9;6;55;mC68;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;S5_{4};Al2H4O9Si2;Nacrite [Al2Si2O5(OH)4, S54]");
    vproto.push_back("A6B36C11_mC212_9_6a_36a_11a;3;106;9;5;163;mC212;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36,x37,y37,z37,x38,y38,z38,x39,y39,z39,x40,y40,z40,x41,y41,z41,x42,y42,z42,x43,y43,z43,x44,y44,z44,x45,y45,z45,x46,y46,z46,x47,y47,z47,x48,y48,z48,x49,y49,z49,x50,y50,z50,x51,y51,z51,x52,y52,z52,x53,y53,z53;-;Cs6O36W11;Cs6W11O36");
    vproto.push_back("ABC_mP6_11_e_e_e;3;6;11;5;10;mP6;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3;-;HO2Y;O(OH)Y");
    vproto.push_back("A3B_mP8_11_3e_e;2;8;11;4;12;mP8;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4;-;Se3Zr;ZrSe3");
    vproto.push_back("A2B5C2_mP18_11_2e_e2f_2e;3;18;11;5;20;mP18;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;K0_{1};K2O5S2;K2S2O5 (K01)");
    vproto.push_back("AB2CD6_mP20_11_e_2e_e_2e2f;4;20;11;6;22;mP20;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7,x8,y8,z8;-;BaC2CaO6;Barytocalcite (BaCa(CO3)2)");
    vproto.push_back("A7B2C2_mP22_11_3e2f_2e_ab;3;22;11;5;20;mP22;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,y8,z8,x9,y9,z9;-;O7Si2Y2;y-Y2Si2O7");
    vproto.push_back("A2B_mC6_12_i_a;2;3;12;4;6;mC6;a,b/a,c/a,beta,x2,z2;-;Cl2Cu;Tolbachite (CuCl2)");
    vproto.push_back("A3B4_mC14_12_ai_2i;2;7;12;4;10;mC14;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4;D7_{a};Ni3Sn4;delta-Ni3Sn4 (D7a)");
    vproto.push_back("AB2C_mC16_12_g_2i_i;3;8;12;5;11;mC16;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,z4;-;FeSe2Tl;Monoclinic FeTlSe2");
    vproto.push_back("ABC3_mC20_12_g_i_ij;3;10;12;5;12;mC20;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,y4,z4;-;MnPS3;MnPS3");
    vproto.push_back("A3BC2_mC24_12_ij_h_gi;3;12;12;5;13;mC24;a,b/a,c/a,beta,y1,y2,x3,z3,x4,z4,x5,y5,z5;B36;H3LiO2;LiOH.H2O (B36)");
    vproto.push_back("AB_mC24_12_3i_3i;2;12;12;4;16;mC24;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;AsSi;SiAs");
    vproto.push_back("A13B4_mC34_12_b6i_2i;2;17;12;4;20;mC34;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9;-;Al13Os4;Os4Al13");
    vproto.push_back("AB6C2D_mC40_12_ad_gh4i_j_bc;4;20;12;6;17;mC40;a,b/a,c/a,beta,y5,y6,x7,z7,x8,z8,x9,z9,x10,z10,x11,y11,z11;-;NiO6Sr2Te;Sr2NiTeO6");
    vproto.push_back("A2B12CD6_mC42_12_i_2i2j_a_ij;4;21;12;6;21;mC42;a,b/a,c/a,beta,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;J1_{7};Cl2H12MgO6;Bischofite (MgCl2.6H2O, J17)");
    vproto.push_back("AB5_mC48_12_2i_ac5i2j;2;24;12;4;24;mC48;a,b/a,c/a,beta,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,y10,z10,x11,y11,z11;D2_{2};MgZn5;D22 (MgZn5?) (problematic)");
    vproto.push_back("AB8C4_mC52_12_i_gi3j_2j;3;26;12;5;24;mC52;a,b/a,c/a,beta,y1,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;S6_{7};AlKO8Si3;Sanidine (KAlSi3O8, S67)");
    vproto.push_back("A2B2C5D24E8_mC82_12_h_i_agh_2i5j_2j;5;41;12;7;34;mC82;a,b/a,c/a,beta,y2,y3,y4,x5,z5,x6,z6,x7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;S4_{2};Ca2H2Mg5O24Si8;Tremolite (Ca2Mg5Si8O22(OH)2, S42)");
    vproto.push_back("A5B2C10D2E2_mC84_12_acghj_bdi_5j_2i_j;5;42;12;7;33;mC84;a,b/a,c/a,beta,y5,y6,x7,z7,x8,z8,x9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;Al5Fe2O10(OH)2Si2;Staurolite (Al5Fe2O10(OH)2Si2)");
    vproto.push_back("A13B4_mC102_12_dg8i5j_4ij;2;51;12;4;47;mC102;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20;-;Al13Fe4;Al13Fe4");
    vproto.push_back("A8B2CD15E2_mC112_12_2i3j_j_ad_g4i5j_2i;5;56;12;7;48;mC112;a,b/a,c/a,beta,y3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20;H4_{23};H8K2MnO12S2;Manganese-leonite [K2Mn(SO4)2.4H2O, H423]");
    vproto.push_back("AB6C11D6E4_mC112_12_e_gi2j_i5j_2i2j_2j;5;56;12;7;46;mC112;a,b/a,c/a,beta,y2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;S4_{5};(H2O)(OH)6O11Mg6Si4;Chrysotile (H4Mg3Si2O9, S45)");
    vproto.push_back("A8B23_mP62_13_4g_c11g;2;62;13;4;49;mP62;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;Mo8O23;High-Temperature Mo8O23");
    vproto.push_back("A2B3C2_mP14_14_e_ae_e;3;14;14;5;13;mP14;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Cl2Hg3O2;HgCl2.2HgO");
    vproto.push_back("ABCD_mP16_14_e_e_e_e;4;16;14;6;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;ClCuHO;Cu(OH)Cl");
    vproto.push_back("AB_mP16_14_2e_2e;2;16;14;4;16;mP16;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-/-;ClI/AsLi;alpha-ICl/LiAs");
    vproto.push_back("A3BC_mP20_14_3e_e_e;3;20;14;5;19;mP20;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Cl3CuK;Sanguite (KCuCl3)");
    vproto.push_back("AB6C2D_mP20_14_a_3e_e_d;4;20;14;6;16;mP20;a,b/a,c/a,beta,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;MnO6Sr2Te;Sr2MnTeO6");
    vproto.push_back("AB6C3_mP20_14_a_3e_de;3;20;14;5;16;mP20;a,b/a,c/a,beta,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;J2_{6};AlF6Na3;Cryolite (Na3AlF6, J26)");
    vproto.push_back("A2B5C4_mP22_14_e_c2e_2e;3;22;14;5;19;mP22;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Cl2O5Sb4;Sb4O5Cl2");
    vproto.push_back("A4B2C4D_mP22_14_2e_e_2e_a;4;22;14;6;19;mP22;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;C4K2N4Ni;K2Ni(CN)4");
    vproto.push_back("A9B2_mP22_14_a4e_e;2;22;14;4;19;mP22;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;D8_{d};Co2Al9;Co2Al9 (D8d)");
    vproto.push_back("A4BC_mP24_14_4e_e_e;3;24;14;5;22;mP24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;O7Si2Y2;gamma-Y2Si2O7");
    vproto.push_back("AB4C_mP24_14_ab_4e_e;3;24;14;5;19;mP24;a,b/a,c/a,beta,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;AuBr4K;Anhydrous KAuBr4");
    vproto.push_back("ABCD3_mP24_14_e_e_e_3e;4;24;14;6;22;mP24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;G0_{12};CHNaO3;Nahcolite (NaHCO3, G012)");
    vproto.push_back("AB_mP24_14_3e_3e;2;24;14;4;22;mP24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;C6Cl6;epsilon-1,2,3,4,5,6-Hexachlorocyclohexane (C6Cl6)");
    vproto.push_back("AB4C_mP24_14_e_4e_e;3;24;14;5;22;mP24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;K4_{1}/-;(NH4)O4S/LaO4P;Ammonium Persulfate [(NH4)2S2O8, K41]/Monazite (LaPO4)");
    vproto.push_back("ABC4_mP24_14_e_e_4e;3;24;14;5;22;mP24;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;H0_{9};AgMnO4;AgMnO4 (H09)");
    vproto.push_back("A2B4C_mP28_14_abe_4e_e;3;28;14;5;22;mP28;a,b/a,c/a,beta,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Cu2O4Se;Monoclinic Cu2OSeO3");
    vproto.push_back("A4BCD_mP28_14_4e_e_e_e;4;28;14;6;25;mP28;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;H0_{10};Cl4IK.H2O;KICl4.H2O (H010)");
    vproto.push_back("A2B3C2D8_mP30_14_e_ce_e_4e;4;30;14;6;25;mP30;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;G7_{4};C2Cu3H2O8;Azurite [Cu3(CO3)2(OH)2, G74]");
    vproto.push_back("A2B5C_mP32_14_2e_5e_ab;3;32;14;5;25;mP32;a,b/a,c/a,beta,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Ca2O5U;Ca2UO5");
    vproto.push_back("A2B5C_mP32_14_2e_5e_e;3;32;14;5;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Gd2O5Si;Gd2SiO5 (RE2SiO5 X1)");
    vproto.push_back("A3B_mP32_14_6e_2e;2;32;14;4;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;O3W;gamma-WO3");
    vproto.push_back("AB4C2D_mP32_14_e_4e_2e_e;4;32;14;6;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;H4_{19};AuBr4(H2O)2K;KAuBr4.2H2O (H419)");
    vproto.push_back("AB_mP32_14_4e_4e;2;32;14;4;28;mP32;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-/B_{l};AsS/AsS;Pararealgar (AsS)/AsS;Realgar (AsS, Bl)");
    vproto.push_back("A7B2C_mP40_14_7e_2e_e;3;40;14;5;34;mP40;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;K6_{2};F7K2Nb;K2NbF7 (K62)");
    vproto.push_back("A6B4C2D6E2FG6_mP54_14_3e_2e_e_3e_e_a_3e;7;54;14;9;43;mP54;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;C6H4K2N6O2PtS6;K2Pt(SCN)6.2H2O");
    vproto.push_back("A11B3_mP56_14_11e_3e;2;56;14;4;46;mP56;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Cs11O3;Cs11O3");
    vproto.push_back("AB3C_mP60_14_3e_9e_3e;3;60;14;5;49;mP60;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;S3_{3} (II);CaO3Si;Parawollastonite (CaSiO3, S33 (II))");
    vproto.push_back("AB20C2D14E2_mP78_14_a_10e_e_7e_e;5;78;14;7;61;mP78;a,b/a,c/a,beta,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20;H4_{4};CuH20N2O14S2;Tutton salt [Cu(NH4)2(SO4)2.6H2O, H44]");
    vproto.push_back("A8B2CD12E2_mP100_14_8e_2e_ad_12e_2e;5;100;14;7;76;mP100;a,b/a,c/a,beta,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26;-;H8K2MnO12S2;Manganese-leonite 110K [K2Mn(SO4)2.4H2O]");
    vproto.push_back("A_mC4_15_e;1;2;15;3;5;mC4;a,b/a,c/a,beta,y1;-;Ga;beta-Ga (obsolete)");
    vproto.push_back("A2B2C_mC20_15_ad_f_e;3;10;15;5;8;mC20;a,b/a,c/a,beta,y3,x4,y4,z4;-;Ag2O2Pb;Ag2PbO2");
    vproto.push_back("A2B_mC24_15_2f_ce;2;12;15;4;11;mC24;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;-;O2Sb;Clinocervantite (beta-Sb2O4)");
    vproto.push_back("AB5C2_mC32_15_e_e2f_f;3;16;15;5;15;mC32;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;NiSe5Ta2;Ta2NiSe5");
    vproto.push_back("AB5CD_mC32_15_e_e2f_e_b;4;16;15;6;13;mC32;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5,x6,y6,z6;S0_{6};CaO5SiTi;Titanite (CaTiSiO5, S06)");
    vproto.push_back("A2BC4D2_mC36_15_f_e_2f_f;4;18;15;6;17;mC36;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;C2(H2O)O4Rb2;Rb2C2O4.H2O");
    vproto.push_back("A7B2C2_mC44_15_e3f_f_f;3;22;15;5;20;mC44;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;O7V2Zn2;alpha-Zn2V2O7");
    vproto.push_back("AB4C6D_mC48_15_e_2f_3f_e;4;24;15;6;21;mC48;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;H4_{6};CaH4O6S;Gypsum (CaSO4.2H2O, H46)");
    vproto.push_back("AB4C4D4E_mC56_15_e_2f_2f_2f_a;5;28;15;7;23;mC56;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;H4_{22};BaC4(H2O)4N4Ni;BaNi(CN)4.4H2O (H422)");
    vproto.push_back("A5BC2_mC64_15_5f_f_2f;3;32;15;5;28;mC64;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;O5SiY2;Y2SiO5 (RE2SiO5 X2)");
    vproto.push_back("AB5CD2_mC72_15_f_5f_f_2f;4;36;15;6;31;mC72;a,b/a,c/a,beta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;S5_{6};AlO5(OH)Si2;Pyrophyllite [AlSi2O5(OH), S56]");
    vproto.push_back("A2BC10D2E4_mC76_15_f_e_5f_f_2f;5;38;15;7;32;mC76;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;S5_{1};Al3KO10(OH)2Si3;Muscovite (KH2Al3Si3O12, S51)");
    vproto.push_back("A2BCD12E3_mC76_15_f_e_b_6f_ef;5;38;15;7;30;mC76;a,b/a,c/a,beta,y2,y3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Fe2MnNaO12P3;Alluaudite [NaMnFe2(PO4)3]");
    vproto.push_back("A3B5C4D2_mC112_15_a3ef_5f_4f_2f;4;56;15;6;43;mC112;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;S5_{5};Mg3O5(OH)4Si2;Chlorite [Mg3(Mg2Al)(Si3Al)O10(OH)8, S55]");
    vproto.push_back("A2B4C2D17E6_mC124_15_f_2f_f_e8f_3f;5;62;15;7;50;mC124;a,b/a,c/a,beta,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;-;BeHNaO8Si3;Eudidymite (BeHNaO8Si3)");
    vproto.push_back("A2B3C9D3E_mC144_15_2f_bcdef_9f_3f_ae;5;72;15;7;51;mC144;a,b/a,c/a,beta,y5,y6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;(H2O)2Na2O9Si3Zr;Catapleiite (Na2ZrSi3O9.2H2O)");
    vproto.push_back("A3B16C20D3_mC168_15_ef_8f_10f_ef;4;84;15;6;66;mC168;a,b/a,c/a,beta,y1,y2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;H4_{20};Cd3H16O20S3;(CdSO4)3.8H2O (H420)");
    vproto.push_back("A8B2CD12E2_mC200_15_8f_2f_ce_2e11f_2f;5;100;15;7;76;mC200;a,b/a,c/a,beta,y2,y3,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27;-;H8K2MnO12S2;Manganese-leonite 185K [K2Mn(SO4)2.4H2O]");
    vproto.push_back("ABC3_oP40_17_abcd_2e_abcd4e;3;40;17;5;29;oP40;a,b/a,c/a,x1,x2,x3,x4,y5,y6,y7,y8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;NaNbO3;NaNbO3");
    vproto.push_back("A2B_oP12_18_2c_c;2;12;18;4;12;oP12;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;O2Te;gamma-TeO2");
    vproto.push_back("AB12C5D2_oP40_18_a_6c_b2c_c;4;40;18;6;32;oP40;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;ClH12N5Zn2;Diamminetriamidodizinc Chloride ([Zn2(NH3)2(NH2)3]Cl)");
    vproto.push_back("AB_oP16_19_2a_2a;2;16;19;4;15;oP16;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;NaP;NaP");
    vproto.push_back("A2B2C_oP20_19_2a_2a_a;3;20;19;5;18;oP20;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;C31;H2O2Zn;W\"ulfingite (epsilon-Zn(OH)2, C31)");
    vproto.push_back("AB4C_oP24_19_a_4a_a;3;24;19;5;21;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;AlCl4Na;NaAlCl4");
    vproto.push_back("A6BC4D_oP48_19_6a_a_4a_a;4;48;19;6;39;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;H6NO4P;Ferroelectric NH4H2PO4");
    vproto.push_back("AB2C_oP80_19_5a_10a_5a;3;80;19;5;63;oP80;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20;-;CH2O;beta-Arabinose [(CH2O)20]");
    vproto.push_back("A14BC11D_oP108_19_14a_a_11a_a;4;108;19;6;84;oP108;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27;H4_{12};H14NiO11S;Morenosite (NiSO4.7H2O, H412)");
    vproto.push_back("AB3_oC16_20_a_bc;2;8;20;4;8;oC16;a,b/a,c/a,x1,y2,x3,y3,z3;D0_{7};CrO3;D07 (CrO3) (obsolete)");
    vproto.push_back("AB4C_oC24_20_b_2c_a;3;12;20;5;11;oC24;a,b/a,c/a,x1,y2,x3,y3,z3,x4,y4,z4;-;AlO4P;AlPO4 ``low cristobalite type''");
    vproto.push_back("AB5C2_oC32_20_b_a2bc_c;3;16;20;5;13;oC32;a,b/a,c/a,x1,y2,y3,y4,x5,y5,z5,x6,y6,z6;K3_{3};AlF5Tl2;Tl2AlF5 (K33)");
    vproto.push_back("A2B7C2_oF88_22_k_bdefghij_k;3;22;22;5;15;oF88;a,b/a,c/a,x3,y4,z5,z6,y7,x8,x9,y9,z9,x10,y10,z10;-;Cd2O7Re2;Predicted Phase IV Cd2Re2O7");
    vproto.push_back("AB6_oP28_29_a_6a;2;28;29;4;24;oP28;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;HgN6;Mercury (II) Azide [Hg(N3)2]");
    vproto.push_back("ABC30DE20F2_oP220_29_a_a_30a_a_20a_2a;6;220;29;8;168;oP220;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36,x37,y37,z37,x38,y38,z38,x39,y39,z39,x40,y40,z40,x41,y41,z41,x42,y42,z42,x43,y43,z43,x44,y44,z44,x45,y45,z45,x46,y46,z46,x47,y47,z47,x48,y48,z48,x49,y49,z49,x50,y50,z50,x51,y51,z51,x52,y52,z52,x53,y53,z53,x54,y54,z54,x55,y55,z55;-;AlCH30NO20S2;Low-Temperature (NH3CH3)Al(SO4)2.12H2O");
    vproto.push_back("A4B7C_oP24_31_2b_a3b_a;3;24;31;5;22;oP24;a,b/a,c/a,y1,z1,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;B4O7Sr;B4SrO7");
    vproto.push_back("A2B6CD8_oP34_31_2a_2a2b_a_4a2b;4;34;31;6;33;oP34;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;H4_{11};Cl2(H2O)6MgO8;Mg(ClO4)2.6H2O (H411)");
    vproto.push_back("A13B4_oP102_31_17a11b_8a2b;2;102;31;4;92;oP102;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8,y9,z9,y10,z10,y11,z11,y12,z12,y13,z13,y14,z14,y15,z15,y16,z16,y17,z17,y18,z18,y19,z19,y20,z20,y21,z21,y22,z22,y23,z23,y24,z24,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36,x37,y37,z37,x38,y38,z38;-;Al13Co4;Orthorhombic Co4Al13");
    vproto.push_back("A17B47_oP128_32_a8c_a23c;2;128;32;4;98;oP128;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33;-;Mo17O47;Mo17O47");
    vproto.push_back("ABC3_oP20_33_a_a_3a;3;20;33;5;18;oP20;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;ILiO3;gamma-LiIO3");
    vproto.push_back("A2B7C2_oP44_33_2a_7a_2a;3;44;33;5;36;oP44;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;Gd2O7Si2;Possible delta-Gd2Si2O7");
    vproto.push_back("A4BCD6_oP48_33_4a_a_a_6a;4;48;33;6;39;oP48;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;B4CsFO6;CsB4O6F");
    vproto.push_back("A2BC4_oP84_33_6a_3a_12a;3;84;33;5;66;oP84;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21;-;B2CaO4;CaB2O4 (III)");
    vproto.push_back("A2BC2_oP10_34_c_a_c;3;10;34;5;10;oP10;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;-;F2-xHxMnOx;MnF2-x(OH)x");
    vproto.push_back("A2BC2_oC20_36_b_a_b;3;10;36;5;11;oC20;a,b/a,c/a,y1,z1,x2,y2,z2,x3,y3,z3;-;N2OSi2;Si2N2O");
    vproto.push_back("A2BC5_oC32_36_b_a_a2b;3;16;36;5;16;oC32;a,b/a,c/a,y1,z1,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;Bi2GeO5;Bi2GeO5");
    vproto.push_back("A4B7C2D2_oC60_36_2b_a3b_2a_b;4;30;36;6;27;oC60;a,b/a,c/a,y1,z1,y2,z2,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;S4_{6};Be4O7(OH)2Si2;Bertrandite (Be4Si2O7(OH)2, S46)");
    vproto.push_back("A3B2_oC80_36_4a4b_2a3b;2;40;36;4;36;oC80;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;Ni3Si2;Ni3Si2");
    vproto.push_back("ABC3_oC80_36_2ab_2ab_2a5b;3;40;36;5;36;oC80;a,b/a,c/a,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;KNO3;alpha-Potassium Nitrate (KNO3) II");
    vproto.push_back("ABC6D15_oC46_38_b_b_2a2d_2ab4d2e;4;23;38;6;26;oC46;a,b/a,c/a,z1,z2,z3,z4,z5,z6,z7,y8,z8,y9,z9,y10,z10,y11,z11,y12,z12,y13,z13,y14,z14,y15,z15;-;FNaNb6O15;NaNb6O15F");
    vproto.push_back("A2B7C2_oC88_40_abc_2b6c_a3b;3;44;40;5;38;oC88;a,b/a,c/a,z1,z2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;Mo2O7Rb2;Rb2Mo2O7");
    vproto.push_back("A5B8CD12_oC104_41_a2b_4b_a_6b;4;52;41;6;41;oC104;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;K3_{5};B5H8KO12;Santite (KB5O8.4H2O, K35)");
    vproto.push_back("A2B_oF24_43_b_a;2;6;43;4;7;oF24;a,b/a,c/a,z1,x2,y2,z2;-;Cs2Se;Cs2Se");
    vproto.push_back("A3B2_oF40_43_ab_b;2;10;43;4;10;oF40;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3;-;Al3Zr2;Zr2Al3");
    vproto.push_back("A2BC4D_oF64_43_b_a_2b_a;4;16;43;6;14;oF64;a,b/a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;H2KO4P;Archerite (KH2PO4)");
    vproto.push_back("A2B7C2_oF88_43_b_a3b_b;3;22;43;5;19;oF88;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;Cu2O7V2;Blossite (alpha-Cu2V2O7)");
    vproto.push_back("A2B4C2D12E3_oF184_43_b_2b_b_6b_ab;5;46;43;7;37;oF184;a,b/a,c/a,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;S6_{10};Al2H4Na2O12Si3;Natrolite (Na2Al2Si3O10.2H2O, S610)");
    vproto.push_back("ABC2_oI8_44_a_a_c;3;4;44;5;7;oI8;a,b/a,c/a,z1,z2,x3,z3;F5_{5};NaNO2;Ferroelectric NaNO2 (F55)");
    vproto.push_back("ABC2_oI8_44_a_a_d;3;4;44;5;7;oI8;a,b/a,c/a,z1,z2,y3,z3;F5_{12};AgNO2;AgNO2 (F512)");
    vproto.push_back("A2B5CD2_oI40_44_2c_abcde_d_e;4;20;44;6;21;oI40;a,b/a,c/a,z1,z2,x3,z3,x4,z4,x5,z5,y6,z6,y7,z7,x8,y8,z8,x9,y9,z9;S2_{2};H2O5SiZn2;Hemimorphite (Zn4Si2O7(OH)2.H2O, S22)");
    vproto.push_back("AB_oI48_44_6d_ab2cde;2;24;44;4;26;oI48;a,b/a,c/a,z1,z2,x3,z3,x4,z4,y5,z5,y6,z6,y7,z7,y8,z8,y9,z9,y10,z10,y11,z11,x12,y12,z12;B30;MgZn;B30 (MgZn?)");
    vproto.push_back("A2B17C6_oI100_46_ab_b8c_3c;3;50;46;5;41;oI100;a,b/a,c/a,z1,y2,z2,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;Nb2O17Zr6;Nb2Zr6O17");
    vproto.push_back("AB2C_oP8_51_e_be_f;3;8;51;5;6;oP8;a,b/a,c/a,z2,z3,z4;-;Bi2Ni3S2;Parkerite (Ni3Bi2S2)");
    vproto.push_back("ABC6D15_oP46_51_f_d_2e2i_aef4i2j;4;46;51;6;24;oP46;a,b/a,c/a,z3,z4,z5,z6,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14,x15,z15;-;FLiNb6O15;LiNb6O15F");
    vproto.push_back("A3B12CDE6_oP276_52_d4e_18e_ce_de_2d8e;5;276;52;7;104;oP276;a,b/a,c/a,z1,x2,x3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27,x28,y28,z28,x29,y29,z29,x30,y30,z30,x31,y31,z31,x32,y32,z32,x33,y33,z33,x34,y34,z34,x35,y35,z35,x36,y36,z36,x37,y37,z37;-;Cl3H12KMgO6;Carnallite [Mg(H2O)6KCl3]");
    vproto.push_back("A2BC_oP16_53_eh_ab_g;3;16;53;5;7;oP16;a,b/a,c/a,x3,y4,y5,z5;F5_{8};F2H5N;NH4HF2 (F58)");
    vproto.push_back("A2BC4D2_oP18_53_h_a_i_e;4;18;53;6;9;oP18;a,b/a,c/a,x2,y3,z3,x4,y4,z4;C45;Cl2CuH4O2;Eriochalcite (CuCl2.2H2O, C45)");
    vproto.push_back("A2B3C8_oP26_55_h_ag_2g2h;3;26;55;5;15;oP26;a,b/a,c/a,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7;-;Nb2Pd3Se8;Nb2Pd3Se8");
    vproto.push_back("A4BCD2_oP32_55_ghi_f_e_gh;4;32;55;6;16;oP32;a,b/a,c/a,z1,z2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;E3_{4};Cl4(H2O)HgK2;K2HgCl4.H2O (E34)");
    vproto.push_back("AB2C5_oP32_55_g_fh_eghi;3;32;55;5;16;oP32;a,b/a,c/a,z1,z2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;-;HoMn2O5;HoMn2O5");
    vproto.push_back("A8B11_oP38_55_g3h_a3g2h;2;38;55;4;21;oP38;a,b/a,c/a,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10;-;B8Ru11;Ru11B8");
    vproto.push_back("A10B3C4_oP68_55_2e2fgh2i_adef_2e2f;3;68;55;5;23;oP68;a,b/a,c/a,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,x13,y13,x14,y14,x15,y15,z15,x16,y16,z16;-;O10Ru3Sr4;Orthorhombic Sr4Ru3O10");
    vproto.push_back("A2BC4_oP56_56_2e_e_4e;3;56;56;5;24;oP56;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;B2CaO4;Calciborite (CaB2O4 II)");
    vproto.push_back("A3B_oP16_57_a2d_d;2;16;57;4;9;oP16;a,b/a,c/a,x2,y2,x3,y3,x4,y4;D0_{10};O3W;D010 (WO3) (em obsolete)");
    vproto.push_back("A4BC_oP24_57_cde_d_a;3;24;57;5;11;oP24;a,b/a,c/a,x2,x3,y3,x4,y4,x5,y5,z5;-;O4SrU;SrUO4");
    vproto.push_back("ABC3_oP40_57_cd_e_cd2e;3;40;57;5;18;oP40;a,b/a,c/a,x1,x2,x3,y3,x4,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;NaNbO3;Lueshite (NaNbO3)");
    vproto.push_back("AB_oP8_58_g_g;2;8;58;4;7;oP8;a,b/a,c/a,x1,y1,x2,y2;-;InS;InS");
    vproto.push_back("A2B3C6_oP22_58_g_af_gh;3;22;58;5;11;oP22;a,b/a,c/a,z2,x3,y3,x4,y4,x5,y5,z5;-;B2Mg3O6;Kotoite (Mg3(BO3)2)");
    vproto.push_back("A4B3_oP28_58_4g_3g;2;28;58;4;17;oP28;a,b/a,c/a,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7;-;In4Se3;In4Se3");
    vproto.push_back("ABC5D2_oP36_58_g_g_3gh_eg;4;36;58;6;19;oP36;a,b/a,c/a,z1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,z8;H2_{7};AsHO5Zn2;Adamite [Zn2(AsO4)(OH), H27]");
    vproto.push_back("A2B7C24D8_oP82_58_g_ae2f_2g5h_2h;4;82;58;6;33;oP82;a,b/a,c/a,z2,z3,z4,x5,y5,x6,y6,x7,y7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;-;H2Mg7O24Si8;Protoanthophyllite (H2Mg7Si8O24)");
    vproto.push_back("A4B2C3_oP18_59_ef_ab_af;3;18;59;5;12;oP18;a,b/a,c/a,z1,z2,z3,y4,z4,x5,z5,x6,z6;G0_{11};H4N2O3;NH4NO3 IV (G011)");
    vproto.push_back("A2BC4_oP28_60_d_c_2d;3;28;60;5;13;oP28;a,b/a,c/a,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4;E3_{2};B2CaO4;CaB2O4 I (E32)");
    vproto.push_back("AB2C6_oP36_60_c_d_3d;3;36;60;5;16;oP36;a,b/a,c/a,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;E5_{1};FeNb2O6;Columbite (FeNb2O6, E51)");
    vproto.push_back("A5B12_oP68_60_c2d_6d;2;68;60;4;28;oP68;a,b/a,c/a,y1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Cr5O12;Cr5O12");
    vproto.push_back("ABC_oP24_61_c_c_c;3;24;61;5;12;oP24;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3;-;CClO;COCl");
    vproto.push_back("A2B4C_oP28_61_c_2c_a;3;28;61;5;12;oP28;a,b/a,c/a,x2,y2,z2,x3,y3,z3,x4,y4,z4;-;Ca2O4Ru;Ca2RuO4");
    vproto.push_back("AB2CD4_oP64_61_c_2c_c_4c;4;64;61;6;27;oP64;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;G7_{2};BBe2HO4;Hambergite [Be2BO3(OH), G72]");
    vproto.push_back("A7BCD_oP80_61_7c_c_c_c;4;80;61;6;33;oP80;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Cl7OPTi;(TiCl4.POCl3)2");
    vproto.push_back("AB3C_oP80_61_2c_6c_2c;3;80;61;5;33;oP80;a,b/a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;S4_{3};MgO3Si;Enstatite (MgSiO3, S43)");
    vproto.push_back("ABCD_oP16_62_c_c_c_c;4;16;62;6;11;oP16;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4;F5_{14};ClBrI(NH4);NH4ClBrI (F514)");
    vproto.push_back("AB3C_oP20_62_c_3c_c;3;20;62;5;13;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5;E2_{4};CdCl3(NH4);NH4CdCl3 (E24)");
    vproto.push_back("ABC3_oP20_62_c_c_cd;3;20;62;5;12;oP20;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4;-/G0_{10}/G0_{2};KNO3/N(NH4)O3/CaCO3;alpha-Potassium Nitrate (KNO3) I/NH4NO3 III (G010)/Aragonite (CaCO3, G02)");
    vproto.push_back("A2B_oP24_62_4c_2c;2;24;62;4;15;oP24;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6;-;Cs2Sb;Cs2Sb");
    vproto.push_back("AB4C_oP24_62_a_2cd_c;3;24;62;5;12;oP24;a,b/a,c/a,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;CuO4S;Chalcocyanite (CuSO4)");
    vproto.push_back("A2B4C_oP28_62_2c_2cd_c;3;28;62;5;16;oP28;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;H1_{6};K2O4S;Arcanite (K2SO4, H16)");
    vproto.push_back("A5BC_oP28_62_3cd_c_c;3;28;62;5;16;oP28;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6;-;O5SV;VOSO4");
    vproto.push_back("AB4C2_oP28_62_c_4c_2c;3;28;62;5;17;oP28;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;E3_{3};FeS4Sb2;Berthierite (FeSb2S4, E33)");
    vproto.push_back("AB6_oP28_62_c_6c;2;28;62;4;17;oP28;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7;-;CuN6;Copper (II) Azide [Cu(N3)2]");
    vproto.push_back("A2B5C_oP32_62_bc_3cd_c;3;32;62;5;16;oP32;a,b/a,c/a,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;S0_{3};Al2O5Si;Sillimanite (Al2SiO5, S03)");
    vproto.push_back("A3B_oP32_62_ab4c_2c;2;32;62;4;15;oP32;a,b/a,c/a,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;O3W;Original beta-WO3 (em obsolete)");
    vproto.push_back("A4BC2D_oP32_62_2cd_b_2c_a;4;32;62;6;14;oP32;a,b/a,c/a,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;E3_{5};Cl4(H2O)K2Sn;K2SnCl4.H2O (E35)");
    vproto.push_back("A4BC2D_oP32_62_2cd_c_d_c;4;32;62;6;17;oP32;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6;-;Cl4(H2O)K2Sn;K2SnCl4.H2O");
    vproto.push_back("A2B2C4D_oP36_62_d_d_2cd_c;4;36;62;6;18;oP36;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;S0_{5};Al2F2O4Si;Topaz (Al2SiO4F2, S05)");
    vproto.push_back("AB2C3D3_oP36_62_c_ac_cd_cd;4;36;62;6;17;oP36;a,b/a,c/a,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7;-;ClCu2H3O3;Atacamite (Cu2(OH)3Cl)");
    vproto.push_back("AB6C2_oP36_62_c_2c2d_d;3;36;62;5;18;oP36;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CaO6Ta2;Rynersonite (Orthorhombic CaTa2O6)");
    vproto.push_back("A2B3C4D_oP40_62_d_cd_2cd_c;4;40;62;6;20;oP40;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;S0_{7};F2Mg3O4Si;Norbergite [Mg(F,OH)2.Mg2SiO4, S07]");
    vproto.push_back("A2B6C3_oP44_62_2c_2c2d_3c;3;44;62;5;23;oP44;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,y8,z8,x9,y9,z9;K5_{1};K2O6S3;K2S3O6 (K51)");
    vproto.push_back("A7B2C2_oP44_62_3c2d_2c_d;3;44;62;5;22;oP44;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;O7Si2Y2;Possible delta-Y2Si2O7");
    vproto.push_back("A8BCD_oP44_62_4c2d_c_c_c;4;44;62;6;23;oP44;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,y8,z8,x9,y9,z9;-;Cl8OPSb;SbCl5.POCl3");
    vproto.push_back("A2BC8D2_oP52_62_d_c_2c3d_d;4;52;62;6;24;oP52;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;S6_{3};B2CaO8Si2;Danburite (CaB2Si2O8, S63)");
    vproto.push_back("A4B3_oP56_62_8c_6c;2;56;62;4;31;oP56;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,z13,x14,z14;-;Mo4P3;Mo4P3");
    vproto.push_back("A3B15C5D_oP96_62_cd_3c6d_3cd_c;4;96;62;6;43;oP96;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16;J1_{8};Cl3H15N5Rh;RhCl2(NH3)5Cl (J18)");
    vproto.push_back("A8B4C4DE8F2_oP108_62_4c2d_2d_2cd_c_4c2d_d;6;108;62;8;49;oP108;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19;F2_{1};C8H4K4MoN8O2;K4[Mo(CN)8].2H2O (F21)");
    vproto.push_back("A4B3_oP112_62_8c4d_4c4d;2;112;62;4;51;oP112;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11,x12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20;-;P4Se3;P4Se3");
    vproto.push_back("ABCD8E3_oP112_62_d_2c_d_4c6d_3d;5;112;62;7;48;oP112;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17;S4_{7};BeHNaO8Si3;Epididymite (BeHNaO8Si3, S47)");
    vproto.push_back("A2B5C22D2E8_oP156_62_d_c2d_2c10d_2c_4d;5;156;62;7;64;oP156;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22;S4_{4};Fe2Mg5O22(OH)2Si8;Anthophyllite (Mg5Fe2Si8O22(OH)2, S44)");
    vproto.push_back("AB22C23D2E2_oP200_62_c_11d_3c10d_d_d;5;200;62;7;80;oP200;a,b/a,c/a,x1,z1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20,x21,y21,z21,x22,y22,z22,x23,y23,z23,x24,y24,z24,x25,y25,z25,x26,y26,z26,x27,y27,z27;-;CaH22O23P2U2;Autunite Ca[(UO2)(PO4)]2(H2O)11");
    vproto.push_back("AB3_oC16_63_c_3c;2;8;63;4;7;oC16;a,b/a,c/a,y1,y2,y3,y4;D0_{d};AsMn3;Mn3As (D0d)");
    vproto.push_back("AB2C2_oC20_63_c_f_2c;3;10;63;5;8;oC20;a,b/a,c/a,y1,y2,y3,y4,z4;E0_{4};FeHO2;Lepidocrocite (gamma-FeO(OH), E04)");
    vproto.push_back("AB2CD_oC20_63_b_f_c_c;4;10;63;6;7;oC20;a,b/a,c/a,y2,y3,y4,z4;-;CFe2SiTh;ThFe2SiC");
    vproto.push_back("A5B_oC24_63_c2f_c;2;12;63;4;9;oC24;a,b/a,c/a,y1,y2,y3,z3,y4,z4;-;Te5Zr;ZrTe5");
    vproto.push_back("A_oC24_63_3f;1;12;63;3;9;oC24;a,b/a,c/a,y1,z1,y2,z2,y3,z3;-;Si24;Si24 Clathrate");
    vproto.push_back("A5B3_oC32_63_cfg_ce;2;16;63;4;10;oC32;a,b/a,c/a,y1,y2,x3,y4,z4,x5,y5;-;Pd5Pu3;Pd5Pu3");
    vproto.push_back("A2B2C6DE2_oC52_63_g_e_fh_c_f;5;26;63;7;14;oC52;a,b/a,c/a,y1,x2,y3,z3,y4,z4,x5,y5,x6,y6,z6;-;Br2Cu2O6PbSe2;Cu2Pb(SeO3)2Br2");
    vproto.push_back("A4BC12D2_oC76_63_eg_c_f3gh_g;4;38;63;6;20;oC76;a,b/a,c/a,y1,x2,y3,z3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,z9;S0_{4};Al4FeO12Si2;S04 (Staurolite, Fe(OH)2Al4Si2O10) (obsolete)");
    vproto.push_back("A10B3C4_oC68_64_2dfg_ad_2d;3;34;64;5;13;oC68;a,b/a,c/a,x2,x3,x4,x5,x6,y7,z7,x8,y8,z8;-;O10Ru3Sr4;Base-centered orthorhombic Sr4Ru3O10");
    vproto.push_back("A2B2C7_oC88_64_ef_df_3f2g;3;44;64;5;21;oC88;a,b/a,c/a,x1,y2,y3,z3,y4,z4,y5,z5,y6,z6,y7,z7,x8,y8,z8,x9,y9,z9;-;Mo2Na2O7;Na2Mo2O7");
    vproto.push_back("A3B8_oC22_65_ag_bd2gh;2;11;65;4;7;oC22;a,b/a,c/a,x4,x5,x6,x7;-;FNb3O7;Nb3O7F");
    vproto.push_back("A2B8CD2_oC26_65_h_r_a_i;4;13;65;6;8;oC26;a,b/a,c/a,x2,y3,x4,y4,z4;E1_{3};Cl2H6MgN2;Mg(NH3)2Cl2 (E13)");
    vproto.push_back("A2BC2D_oC24_67_m_a_n_g;4;12;67;6;8;oC24;a,b/a,c/a,z2,y3,z3,x4,z4;F5_{7};H2(NH4)O2P;NH4H2PO2 (F57)");
    vproto.push_back("A2B4C_oF56_70_g_h_a;3;14;70;5;7;oF56;a,b/a,c/a,z2,x3,y3,z3;H1_{7};Na2O4S;Thenardite [Na2SO4 (V), H17]");
    vproto.push_back("ABC2_oI16_71_g_i_eh;3;8;71;5;7;oI16;a,b/a,c/a,x1,y2,y3,z4;-;CsFeS2;CsFeS2 (100K)");
    vproto.push_back("AB6C3_oI20_71_a_in_cj;3;10;71;5;7;oI20;a,b/a,c/a,z3,z4,x5,y5;-;AlF6Na3;High-Temperature Cryolite (Na3AlF6)");
    vproto.push_back("A2B5_oI28_72_j_bfj;2;14;72;4;8;oI28;a,b/a,c/a,x2,x3,y3,x4,y4;D8_{g};Ga2Mg5;Ga2Mg5 (D8g)");
    vproto.push_back("ABC4D_oI28_74_a_d_hi_e;4;14;74;6;8;oI28;a,b/a,c/a,z3,y4,z4,x5,z5;-;CuLiO4V;LiCuVO4");
    vproto.push_back("A2B6C2D_oI44_74_h_ij_i_e;4;22;74;6;13;oI44;a,b/a,c/a,z1,y2,z2,x3,z3,x4,z4,x5,y5,z5;E1_{2};Cl2H6N2Zn;Zn(NH3)2Cl2 (E12)");
    vproto.push_back("A4B2C3_tP72_77_8d_ab2c2d_6d;3;72;77;5;54;tP72;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19,x20,y20,z20;-;H4N2O3;Gwihabaite [NH4NO3 (V)]");
    vproto.push_back("AB5C_tP14_85_c_cg_b;3;14;85;5;7;tP14;a,c/a,z2,z3,x4,y4,z4;-;MoO5P;MoPO5");
    vproto.push_back("A3B6CD_tP44_85_bcg_3g_ac_e;4;44;85;6;16;tP44;a,c/a,z3,z4,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;E2_{6};(Br,Cl)3(H2O)6KMg;Bromocarnallite (KMg(H2O)6(Cl,Br)3, E26)");
    vproto.push_back("A2BC_tP32_86_2g_g_g;3;32;86;5;14;tP32;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;E1_{4};Cl2NP;PNCl2 (E14)");
    vproto.push_back("AB6C_tP32_86_d_3g_c;3;32;86;5;11;tP32;a,c/a,x3,y3,z3,x4,y4,z4,x5,y5,z5;J1_{11};Na(OH)6Sb;NaSb(OH)6 (J111)");
    vproto.push_back("ABC3_tP40_86_g_g_3g;3;40;86;5;17;tP40;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;ILiO3;beta-LiIO3");
    vproto.push_back("A4B11C2_tP68_86_2g_ab5g_g;3;68;86;5;26;tP68;a,c/a,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Nd4O11Re2;Nd4Re2O11");
    vproto.push_back("AB6C2D_tI20_87_a_eh_d_b;4;10;87;6;5;tI20;a,c/a,z4,x5,y5;-;NiO6Sr2W;Sr2NiWO6");
    vproto.push_back("AB4C24D12_tI82_87_a_h_2h2i_hi;4;41;87;6;19;tI82;a,c/a,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,z6,x7,y7,z7,x8,y8,z8;S6_{4};ClNa4O24(Al3Si9);Marialite Scapolite [Na4Cl(AlSi3)3O24, S64]");
    vproto.push_back("A9B4C20_tI132_88_a2f_f_5f;3;66;88;5;26;tI132;a,c/a,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9;-;Ge9Na4O20;Na4Ge9O20");
    vproto.push_back("AB2C3D2_tP16_90_c_f_ce_e;4;16;90;6;7;tP16;a,c/a,z1,z2,x3,x4,x5;G7_{5};CCl2O3Pb2;G75 (PbCO3.PbCl2, Phosgenite) (obsolete)");
    vproto.push_back("A12BC10D_tP96_92_6b_a_5b_a;4;96;92;6;37;tP96;a,c/a,x1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;H4_{5};H12NiO10S;Retgersite (alpha-NiSO4.6H2O, H45)");
    vproto.push_back("A2B7C2_tI44_98_f_bcde_f;3;22;98;5;7;tI44;a,c/a,z2,x3,x4,x5,x6;-;Cd2O7Re2;Phase III Cd2Re2O7");
    vproto.push_back("ABC2_tP8_100_b_a_c;3;8;100;5;6;tP8;a,c/a,z1,z2,x3,z3;F5_{4};Cl(NH4)O2;F54 (NH4ClO2) (obsolete)");
    vproto.push_back("ABC3_tP10_100_b_a_bc;3;10;100;5;7;tP10;a,c/a,z1,z2,z3,x4,z4;G0_{9};N(NH4)O3;NH4NO3 II (G09)");
    vproto.push_back("A6B2C_tP72_103_abc5d_2d_abc;3;72;103;5;29;tP72;a,c/a,z1,z2,z3,z4,z5,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;-;O6Se2V;VSe2O6");
    vproto.push_back("A2BC4_tP7_111_f_a_n;3;7;111;5;4;tP7;a,c/a,x3,z3;E3_{1};Ag2HgI4;E31 (beta-Ag2HgI4) (obsolete)");
    vproto.push_back("AB4CD2_tP16_113_c_f_a_e;4;16;113;6;8;tP16;a,c/a,z2,x3,z3,x4,y4,z4;-;ClH4NO2;Ammonium Chlorite (NH4ClO2)");
    vproto.push_back("A2B12C4D4E_tP46_114_d_3e_e_e_a;5;46;114;7;18;tP46;a,c/a,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;H4_{17};Ag2H12N4O4S;Ag2SO4.4NH3 (H417)");
    vproto.push_back("A19B15_tP68_114_bc4e_ac3e;2;68;114;4;25;tP68;a,c/a,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;C19Sc15;C19Sc15");
    vproto.push_back("ABC2_tP4_115_a_c_g;3;4;115;5;3;tP4;a,c/a,z3;F6_{1};CuFeS2;F61 (Chalcopyrite, CuFeS2) (obsolete)");
    vproto.push_back("AB2C_tI8_119_c_e_a;3;4;119;5;3;tI8;a,c/a,z3;-;FeS2Tl;Tetragonal TlFeS2");
    vproto.push_back("A2B7C2_tI44_119_i_bdefgh_i;3;22;119;5;10;tI44;a,c/a,z3,z4,x5,x6,x7,z7,x8,z8;-;Cd2O7Re2;Phase II Cd2Re2O7");
    vproto.push_back("AB8C8D_tI72_120_c_2i_2i_b;4;36;120;6;14;tI72;a,c/a,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;H4_{3};BeH8O8S;BeSO4.4H2O (H43)");
    vproto.push_back("AB3_tI32_121_g_f2i;2;16;121;4;8;tI32;a,c/a,x1,x2,x3,z3,x4,z4;-;SV3;alpha-V3S");
    vproto.push_back("A2B2C6D_tI44_121_i_i_ij_c;4;22;121;6;11;tI44;a,c/a,x2,z2,x3,z3,x4,z4,x5,y5,z5;-;B2Cu2O6Sr;SrCu2(BO3)2");
    vproto.push_back("A2BC2_tI40_122_e_d_e;3;20;122;5;9;tI40;a,c/a,x1,x2,y2,z2,x3,y3,z3;F1_{1};C2HgN2;Mercury Cyanide [Hg(CN)2, F11]");
    vproto.push_back("A4BC4D_tI40_122_e_b_e_a;4;20;122;6;8;tI40;a,c/a,x3,y3,z3,x4,y4,z4;H2_{2};H2KO4P;KH2PO4 (H22)");
    vproto.push_back("AB2_tI48_122_cd_2e;2;24;122;4;10;tI48;a,c/a,z1,x2,x3,y3,z3,x4,y4,z4;-;NaS2;NaS2");
    vproto.push_back("A8BC4D_tI56_122_2e_b_e_a;4;28;122;6;11;tI56;a,c/a,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;H6NO4P;NH4H2PO4");
    vproto.push_back("A3BC_tP5_123_cg_a_d;3;5;123;5;3;tP5;a,c/a,z4;E2_{5};Cl3Hg(NH4);NH4HgCl3 (E25)");
    vproto.push_back("AB4C_tP6_123_d_eh_a;3;6;123;5;3;tP6;a,c/a,z4;H0_{8};AlF4Tl;TlAlF4 (H08)");
    vproto.push_back("A4B2C_tP7_123_j_e_a;3;7;123;5;3;tP7;a,c/a,x3;H1_{5};Cl4K2Pt;K2PtCl4 (H15)");
    vproto.push_back("A8B2C_tP11_123_r_f_a;3;11;123;5;4;tP11;a,c/a,x3,z3;E6_{1};Sr(OH)2(H2O)8;E61 (Sr(OH)2(H2O)8) (obsolete)");
    vproto.push_back("A8B2C_tP11_123_r_h_a;3;11;123;5;5;tP11;a,c/a,z2,x3,z3;E6_{2};(H2O)8SrO2;E62 [SrO2(H2O)8] (empossibly obsolete)");
    vproto.push_back("AB8C2_tP22_124_a_n_h;3;22;124;5;6;tP22;a,c/a,z2,x3,y3,z3;-;Ca(H2O)8O2;CaO2(H2O)8");
    vproto.push_back("ABC4D2E8_tP32_126_a_b_h_e_k;5;32;126;7;7;tP32;a,c/a,z3,x4,x5,y5,z5;J1_{9};AgCoN4(NH3)2O8;Ag[Co(NH3)2(NO2)4] (J19)");
    vproto.push_back("A4B10C2D34E4F9_tP252_126_k_ce2k_f_h8k_k_d2k;6;252;126;8;46;tP252;a,c/a,z3,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15,x16,y16,z16,x17,y17,z17,x18,y18,z18,x19,y19,z19;S2_{3};Al4Ca10(Mg,Fe)2O34(OH)4Si9;Vesuvianite (Ca10Al4(Mg,Fe)2Si9O34(OH)4, S23)");
    vproto.push_back("A2BC4D_tP16_127_h_d_i_a;4;16;127;6;5;tP16;a,c/a,x3,x4,y4;H4_{9};Cl2(H2O)(NH3)4Pd;Pd(NH3)4Cl2.H2O (H49)");
    vproto.push_back("AB2C3D2_tP32_127_g_eh_gk_k;4;32;127;6;10;tP32;a,c/a,z1,x2,x3,x4,x5,z5,x6,z6;-;CCl2O3Pb2;Phosgenite [Pb2Cl2(CO3)]");
    vproto.push_back("A3B14C5_tP44_128_ac_ehi_bg;3;44;128;5;9;tP44;a,c/a,z4,x5,x6,y6,x7,y7,z7;K7_{5};Al3F14Na5;Chiolite (Na5Al3F14, K75)");
    vproto.push_back("A4BC16DE28F8_tP116_128_h_a_2i_b_g3i_i;6;116;128;8;23;tP116;a,c/a,x3,x4,y4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;S5_{2};Ca4FH16KO28Si8;Apophyllite (KCa4Si8O20F.8H2O, S52)");
    vproto.push_back("ABCD_tP8_129_b_c_a_c;4;8;129;6;4;tP8;a,c/a,z3,z4;-;AgLaOS;LaOAgS");
    vproto.push_back("AB4C_tP12_129_c_i_a;3;12;129;5;5;tP12;a,c/a,z2,y3,z3;B25;BrH4N;NH4Br (B25)");
    vproto.push_back("AB4C6DE_tP26_129_c_j_2ci_a_c;5;26;129;7;10;tP26;a,c/a,z2,z3,z4,z5,y6,z6,x7,z7;H5_{10};Ca(H2O)6O12P2U2;Meta-autunite (I) [Ca(UO2)2(PO4)2.6H2O, H510]");
    vproto.push_back("A18B10C_tP116_130_2c4g_2c2g_a;3;116;130;5;24;tP116;a,c/a,z2,z3,z4,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;-;H18O10Sr;Sr(OH)2(H2O)8");
    vproto.push_back("A4B2C_tP14_136_i_g_b;3;14;136;5;5;tP14;a,c/a,x2,x3,y3;-;Fe4Si2Zr;ZrFe4Si2");
    vproto.push_back("A2B3_tP20_136_j_dfg;2;20;136;4;6;tP20;a,c/a,x2,x3,x4,z4;-;Al2Zr3;Zr3Al2");
    vproto.push_back("A4BC4D2E2_tP26_136_fg_a_j_d_e;5;26;136;7;7;tP26;a,c/a,z3,x4,x5,x6,z6;H4_{1};Cl4CuH4K2O2;K2CuCl4.2H2O (H41)");
    vproto.push_back("AB14C2_tP68_136_f_ce2j2k_fg;3;68;136;5;16;tP68;a,c/a,z2,x3,x4,x5,x6,z6,x7,z7,x8,y8,z8,x9,y9,z9;-;BFe14Nd2;Nd2Fe14B");
    vproto.push_back("A2B3_tI10_139_e_ae;2;5;139;4;4;tI10;a,c/a,z2,z3;-;Au2Nb3;Au2Nb3");
    vproto.push_back("A4B2C2D_tI18_139_h_d_e_a;4;9;139;6;4;tI18;a,c/a,z3,x4;J1_{5};Cl4K2O2Os;K2OsO2Cl4 (J15)");
    vproto.push_back("AB3C_tI20_139_ab_eh_d;3;10;139;5;4;tI20;a,c/a,z4,x5;K7_{6};AuCl3Cs;AuCsCl3 (K76)");
    vproto.push_back("A10B4C3_tI34_139_c2eg_2e_ae;3;17;139;5;8;tI34;a,c/a,z3,z4,z5,z6,z7,z8;-;O10Sr4Ti3;Sr4Ti3O10");
    vproto.push_back("A6B2C3D_tI168_139_egikl2m_ejn_bh2n_acf;4;84;139;6;21;tI168;a,c/a,z4,z5,z7,x8,x9,x10,x11,x12,y12,x13,z13,x14,z14,y15,z15,y16,z16,y17,z17;J3_{1};Cl6(H2O)2K3Tl;K3TlCl6.2H2O (J31)");
    vproto.push_back("A2BC_tI16_140_h_d_a;3;8;140;5;3;tI16;a,c/a,x3;F5_{2};F2HK;KHF2 (F52)");
    vproto.push_back("AB6_tI28_140_a_hk;2;14;140;4;5;tI28;a,c/a,x2,x3,y3;D2_{c};MnU6;U6Mn (D2c)");
    vproto.push_back("A5BC2_tI32_140_bl_a_h;3;16;140;5;5;tI32;a,c/a,x3,x4,z4;K3_{4};Br5(NH4)Pb2;NH4Pb2Br5 (K34)");
    vproto.push_back("A5BC3_tI36_140_cl_b_ah;3;18;140;5;5;tI36;a,c/a,x4,x5,z5;K3_{1};Cl5CoCs3;Cs3CoCl5 (K31)");
    vproto.push_back("A31B20_tI204_140_b2gh3m_ac2fh3l;2;102;140;4;23;tI204;a,c/a,z4,z5,z6,z7,x8,x9,x10,z10,x11,z11,x12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;Pu31Rh20;Pu31Rh20");
    vproto.push_back("AB11_tI48_141_a_bdi;2;24;141;4;5;tI48;a,c/a,x4,y4,z4;-;BaCd11;BaCd11");
    vproto.push_back("A2B3_tI160_142_deg_3g;2;80;142;4;16;tI160;a,c/a,z1,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;As2Cd3;Cd3As2");
    vproto.push_back("A2B2C3D12E4_tI184_142_f_f_be_3g_g;5;92;142;7;17;tI184;a,c/a,x2,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8;S6_{1};Al(H2O)NaO6Si2;Analcime (NaAlSi2O6.H2O, S61)");
    vproto.push_back("AB3C9D_hP28_143_2a_2d_6d_bc;4;28;143;6;30;hP28;a,c/a,z1,z2,z3,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12;-;BLa3O9W;La3BWO9 (P3)");
    vproto.push_back("AB3C_hP45_144_3a_9a_3a;3;45;144;5;47;hP45;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14,x15,y15,z15;-;NO3Rb;RbNO3 (IV)");
    vproto.push_back("A2B3C_hP12_147_abd_g_d;3;12;147;5;7;hP12;a,c/a,z3,z4,x5,y5,z5;G3_{2};Na2O3S;Na2SO3 (G32)");
    vproto.push_back("A2BCD6_hR10_148_c_a_b_f;4;10;148;6;6;hR10;a,c/a,x3,x4,y4,z4;G1_{1};C2CaMgO6;Dolomite [MgCa(CO3)2, G11]");
    vproto.push_back("A6B6CD_hR14_148_f_f_b_a;4;14;148;6;8;hR14;a,c/a,x3,y3,z3,x4,y4,z4;I6_{1};Cl6(H2O)6NiSn;Ni(H2O)6SnCl6 (I61)");
    vproto.push_back("A6B2C6D_hR15_148_f_c_f_a;4;15;148;6;9;hR15;a,c/a,x2,x3,y3,z3,x4,y4,z4;H6_{2};H6K2O6Sn;K2Sn(OH)6 (H62)");
    vproto.push_back("ABC8D2_hP12_150_b_a_dg_d;4;12;150;6;7;hP12;a,c/a,z3,z4,x5,y5,z5;H3_{2};AlKO8S2;Steklite [KAl(SO4)2, H32]");
    vproto.push_back("A2B9C3_hP14_150_d_eg_ad;3;14;150;5;8;hP14;a,c/a,z2,z3,x4,x5,y5,z5;K7_{3};As2Cl9Cs3;Cs3As2Cl9 (K73)");
    vproto.push_back("A2B12C6D_hP21_150_d_2g_ef_a;4;21;150;6;11;hP21;a,c/a,z2,x3,x4,x5,y5,z5,x6,y6,z6;-;Cl2H12O6Sr;SrCl2.(H2O)6");
    vproto.push_back("AB2CD6_hP30_150_e_c2d_f_3g;4;30;150;6;16;hP30;a,c/a,z1,z2,z3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;BaC2CaO6;Paralstonite (BaCa(CO3)2)");
    vproto.push_back("AB3C_hP30_150_ef_3g_c2d;3;30;150;5;16;hP30;a,c/a,z1,z2,z3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8;K1_{1};KO3S;KSO3 (K11)");
    vproto.push_back("AB2C2DE3_hR9_155_b_c_c_a_e;5;9;155;7;5;hR9;a,c/a,x3,x4,y5;-;BBe2F2KO3;KBe2BO3F2");
    vproto.push_back("ABC3_hR5_160_a_a_b;3;5;160;5;6;hR5;a,c/a,x1,x2,x3,z3;G0_{7}/-;BrKO3/KNO3;KBrO3 (G07)/gamma-Potassium Nitrate (KNO3)");
    vproto.push_back("AB3C2D_hR7_160_a_b_2a_a;4;7;160;6;8;hR7;a,c/a,x1,x2,x3,x4,x5,z5;S5_{7};FeO3(OH)2Si;Cronstedtite Fe(Fe,Si)[(OH)2,O]O3, S57");
    vproto.push_back("A3B7C_hR11_160_b_a2b_a;3;11;160;5;10;hR11;a,c/a,x1,x2,x3,z3,x4,z4,x5,z5;-;Fe3O7P;Fe3PO7");
    vproto.push_back("AB4C8_hR13_160_a_ab_2a2b;3;13;160;5;12;hR13;a,c/a,x1,x2,x3,x4,x5,z5,x6,z6,x7,z7;-;GaMo4S8;Low-Temperature GaMo4S8");
    vproto.push_back("A3B24C_hR28_160_b_2b3c_a;3;28;160;5;18;hR28;a,c/a,x1,x2,z2,x3,z3,x4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;I3S24Sb;SbI3S24");
    vproto.push_back("A2BC4_hR42_161_2b_b_4b;3;42;161;5;23;hR42;a,c/a,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;B2BaO4;alpha-BaB2O4 (Low-Temperature)");
    vproto.push_back("A2B6C_hP9_162_d_k_a;3;9;162;5;4;hP9;a,c/a,x3,z3;I1_{3};Cl2(H2O)6Sr;I13 (SrCl2.(H2O)6) (obsolete)");
    vproto.push_back("A6BC2_hP9_162_k_a_d;3;9;162;5;4;hP9;a,c/a,x3,z3;-;O6PbSb2;Rosiaite (PbSb2O6)");
    vproto.push_back("A6BC_hP16_163_i_b_c;3;16;163;5;5;hP16;a,c/a,x3,y3,z3;J1_{12};F4NaSb;NaSbF4(OH)2 (J112)");
    vproto.push_back("ABC6D_hP18_163_d_b_i_c;4;18;163;6;5;hP18;a,c/a,x4,y4,z4;-;AlCaF6Li;Colquiriite (LiCaAlF6)");
    vproto.push_back("AB3_hP4_164_b_ad;2;4;164;4;3;hP4;a,c/a,z3;D0_{13};AlCl3;D013 (AlCl3) (obsolete)");
    vproto.push_back("A2BC2_hP5_164_d_a_d;3;5;164;5;4;hP5;a,c/a,z2,z3;-;H2MgO2;Brucite [Mg(OH)2]");
    vproto.push_back("A2BC6_hP9_164_d_a_i;3;9;164;5;5;hP9;a,c/a,z2,x3,z3;H6_{3};K2Pt(SCN)6;K2Pt(SCN)6 (H63)");
    vproto.push_back("A6B2C_hP9_164_i_d_a;3;9;164;5;5;hP9;a,c/a,z2,x3,z3;J1_{6};F6(NH4)2Si;Bararite (Trigonal (NH4)2SiF6, J16)");
    vproto.push_back("A6BC2_hP9_164_i_a_d;3;9;164;5;5;hP9;a,c/a,z2,x3,z3;J1_{13};F6GeK2;K2GeF6 (J113)");
    vproto.push_back("AB2C3_hP12_164_d_ae_i;3;12;164;5;5;hP12;a,c/a,z2,x4,z4;-;HgPt2Se3;Jacutingaite (Pt2HgSe3)");
    vproto.push_back("AB_hP12_164_c2d_c2d;2;12;164;4;8;hP12;a,c/a,z1,z2,z3,z4,z5,z6;-;BiSe;Nevskite (BiSe)");
    vproto.push_back("A16B2C_hP19_164_2d2i_d_b;3;19;164;5;9;hP19;a,c/a,z2,z3,z4,x5,z5,x6,z6;-;H16Li2Mg;Predicted Li2MgH16 300GPa");
    vproto.push_back("AB4C_hR6_166_b_2c_a;3;6;166;5;4;hR6;a,c/a,x3,x4;-;CaO4U;CaUO4");
    vproto.push_back("AB4C2_hR7_166_a_2c_c;3;7;166;5;5;hR7;a,c/a,x2,x3,x4;-;CaCu4P2;CaCu4P2");
    vproto.push_back("ABC6_hR8_166_a_b_h;3;8;166;5;4;hR8;a,c/a,x3,z3;-;KNO3;beta-Potassium Nitrate (KNO3)");
    vproto.push_back("A13B2_hR15_166_b2h_c;2;15;166;4;7;hR15;a,c/a,x2,x3,z3,x4,z4;D1_{g};B13C2;B13C2 ``B4C'' (D1g)");
    vproto.push_back("AB4C2_hR28_166_2c_2c2h_abh;3;28;166;5;12;hR28;a,c/a,x3,x4,x5,x6,x7,z7,x8,z8,x9,z9;-;CuS4Ti2;Rhombohedral CuTi2S4");
    vproto.push_back("A5B21C24D12_hR62_166_a2c_ehi_fg2h_i;4;62;166;6;18;hR62;a,c/a,x2,x3,x5,x6,x7,z7,x8,z8,x9,z9,x10,y10,z10,x11,y11,z11;S3_{4} (I);(Ca1.4,Sr0.3)(H2O)13O24(Si8.3,Al3.8);Chabazite (Ca1.4Sr0.3Al3.8Si8.3O24.13H2O, S34 (I))");
    vproto.push_back("A3B_hR8_167_e_b;2;8;167;4;3;hR8;a,c/a,x2;D0_{12};F3Fe;FeF3 (D012)");
    vproto.push_back("A3BC6_hR20_167_e_b_f;3;20;167;5;6;hR20;a,c/a,x2,x3,y3,z3;J2_{2};Cl3Cr(H2O)6;CrCl3(H2O)6 (J22)");
    vproto.push_back("A6BC3D_hR22_167_f_b_e_a;4;22;167;6;6;hR22;a,c/a,x3,x4,y4,z4;-;Cl6FeK3Na;Rinneite (K3NaFeCl6)");
    vproto.push_back("A9B3C2_hR28_167_ef_e_c;3;28;167;5;8;hR28;a,c/a,x1,x2,x3,x4,y4,z4;K7_{2};Cl9Cs3Tl2;Cs3Tl2Cl9 (K72)");
    vproto.push_back("A2BC4_hR42_167_f_ac_2f;3;42;167;5;12;hR42;a,c/a,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;B2BaO4;beta-BaB2O4 (High-Temperature)");
    vproto.push_back("A25B21_hR92_167_b2e3f_e3f;2;92;167;4;23;hR92;a,c/a,x2,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;-;Re25Zr21;Zr21Re25");
    vproto.push_back("ABC3_hP10_173_b_a_c;3;10;173;5;7;hP10;a,c/a,z1,z2,x3,y3,z3;-;ILiO3;alpha-LiIO3");
    vproto.push_back("ABC4D_hP14_173_a_b_bc_b;4;14;173;6;9;hP14;a,c/a,z1,z2,z3,z4,x5,y5,z5;H1_{4};KLiO4S;LiKSO4 (H14)");
    vproto.push_back("AB3C7D_hP24_173_a_c_b2c_b;4;24;173;6;14;hP24;a,c/a,z1,z2,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;CuLa3S7Si;La3CuSiS7");
    vproto.push_back("AB3C9D_hP28_173_a_c_3c_b;4;28;173;6;16;hP28;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6;-;BLa3O9W;La3BWO9 (P63)");
    vproto.push_back("A3BCD3E15F3_hP52_173_c_b_b_c_5c_c;6;52;173;8;28;hP52;a,c/a,z1,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10;S3_{3} (I);Al3CCaNa3O15Si3;Crancrinite (Na6Ca2Al6Si6O24(CO3)2, S33 (I))");
    vproto.push_back("A3B2_hP20_176_2h_ah;2;20;176;4;8;hP20;a,c/a,x2,y2,x3,y3,x4,y4;D8_{k};S12Th7;Th7S12 (D8k)");
    vproto.push_back("A9B3C2_hP28_176_hi_af_f;3;28;176;5;9;hP28;a,c/a,z2,z3,x4,y4,x5,y5,z5;K7_{1};Cl9K3W2;K3W2Cl9 (K71)");
    vproto.push_back("A10B7_hP34_176_c3h_b2h;2;34;176;4;12;hP34;a,c/a,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7;-;Rh20Si13;Rh20Si13");
    vproto.push_back("A9B2C9_hP40_176_hi_f_hi;3;40;176;5;13;hP40;a,c/a,z1,x2,y2,x3,y3,x4,y4,z4,x5,y5,z5;F4_{1};C9Fe2O9;Fe2(CO)9 (F41)");
    vproto.push_back("A5BC12D3_hP42_176_fh_a_2hi_h;4;42;176;6;14;hP42;a,c/a,z2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,z7;H5_{7};Ca5FO12P3;Fluorapatite [Ca5F(PO4)3, H57]");
    vproto.push_back("A2BCD2_hP18_180_f_c_b_i;4;18;180;6;4;hP18;a,c/a,z3,x4;-;Hg2INaO2;Hg2O2NaI");
    vproto.push_back("ABC3_hP10_182_c_b_g;3;10;182;5;3;hP10;a,c/a,x3;E2_{3};ILiO3;E23 (LiIO3) (obsolete)");
    vproto.push_back("A2BC6_hP18_182_f_b_gh;3;18;182;5;5;hP18;a,c/a,z2,x3,x4;H2_{8};Al2BaO4;BaAl2O4 (H28)");
    vproto.push_back("AB2_hP6_186_b_ab;2;6;186;4;5;hP6;a,c/a,z1,z2,z3;C27;CdI2;C27 (CdI2) (emquestionable)");
    vproto.push_back("ABCD_hP8_186_b_b_a_a;4;8;186;6;6;hP8;a,c/a,z1,z2,z3,z4;E0_{3};CdClHO;Cd(OH)Cl (E03)");
    vproto.push_back("A3B8C2_hP26_186_c_ab2c_2b;3;26;186;5;12;hP26;a,c/a,z1,z2,z3,z4,x5,z5,x6,z6,x7,z7;-;Mo3O8Zn2;Zn2Mo3O8");
    vproto.push_back("A4BC7D_hP26_186_ac_b_a2c_b;4;26;186;6;12;hP26;a,c/a,z1,z2,z3,z4,x5,z5,x6,z6,x7,z7;E9_{2};Be4NaO7Sb;Swedenborgite (NaBe4SbO7, E92)");
    vproto.push_back("AB6CD7_hP30_186_b_d_a_b2c;4;30;186;6;12;hP30;a,c/a,z1,z2,z3,x4,z4,x5,z5,x6,y6,z6;H4_{18};ClH6LiO7;LiClO4.3H2O (H418)");
    vproto.push_back("A3B9CD9_hP44_186_c_3c_b_cd;4;44;186;6;16;hP44;a,c/a,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,y7,z7;G2_{2};Br3(H2O)9NdO9;Nd(BrO3)3.9H2O (G22)");
    vproto.push_back("A3B3C2_hP16_187_jk_jk_ck;3;16;187;5;7;hP16;a,c/a,x2,x3,x4,x5,x6;-;As3Cr3K2;Cr-233 Quasi-One-Dimensional Superconductor (K2Cr3As3)");
    vproto.push_back("A7B_hP24_187_ai2j2kn_j;2;24;187;4;10;hP24;a,c/a,z2,x3,x4,x5,x6,x7,x8,z8;-;Cs7O;Cs7O");
    vproto.push_back("AB3C_hP20_190_ac_i_f;3;20;190;5;6;hP20;a,c/a,z3,x4,y4,z4;K1_{2};CsO3S;CsSO3 (K12)");
    vproto.push_back("ABCD3_hP36_190_h_g_af_hi;4;36;190;6;11;hP36;a,c/a,z2,x3,x4,y4,x5,y5,x6,y6,z6;G7_{1};CCeFO3;Bastn\"asite [CeF(CO3)]");
    vproto.push_back("A12B_hP13_191_cdei_a;2;13;191;4;4;hP13;a,c/a,z4,z5;D2_{a};Be12Ti;TiBe12 (approximate, D2a)");
    vproto.push_back("A4B5_hP18_193_bg_dg;2;18;193;4;4;hP18;a,c/a,x3,x4;-;Ga4Ti5;Ti5Ga4");
    vproto.push_back("A3B_hP24_193_ack_g;2;24;193;4;5;hP24;a,c/a,x3,x4,z4;D0_{6};F3La;D06 (Tysonite, LaF3) (obsolete)");
    vproto.push_back("A2B3_hP10_194_f_bf;2;10;194;4;4;hP10;a,c/a,z2,z3;D5_{b};Pt2Sn3;Pt2Sn3 (D5b)");
    vproto.push_back("AB2C2_hP10_194_a_bc_f;3;10;194;5;3;hP10;a,c/a,z4;-;CoNa0.74O2;Na0.74CoO2");
    vproto.push_back("AB2C2_hP10_194_a_f_f;3;10;194;5;4;hP10;a,c/a,z2,z3;-;EuIn2P2;EuIn2P2");
    vproto.push_back("A10B_hP22_194_bhj_c;2;22;194;4;5;hP22;a,c/a,x3,x4,y4;-;H10Hf;Proposed 300GPa HfH10");
    vproto.push_back("AB3C2_hP24_194_f_k_bh;3;24;194;5;6;hP24;a,c/a,z2,x3,x4,z4;-;CoGa3Lu2;Lu2CoGa3");
    vproto.push_back("A9B2C3_hP28_194_hk_f_bf;3;28;194;5;7;hP28;a,c/a,z2,z3,x4,x5,z5;-;Cl9Cr2Cs3;Cs3Cr2Cl9");
    vproto.push_back("A3B2C9D3E_hP36_194_g_f_hk_h_a;5;36;194;7;7;hP36;a,c/a,z2,x4,x5,x6,z6;S3_{4} (II);(H2O)Na2O9Si3Zr;S34 (II) (Catapleiite, Na2Zr(SiO3)3.H2O) (obsolete)");
    vproto.push_back("A2B3_hP60_194_3fk_cdef2k;2;60;194;4;13;hP60;a,c/a,z3,z4,z5,z6,z7,x8,z8,x9,z9,x10,z10;D5_{6};Al2O3;beta-Alumina (Al2O3, D56)");
    vproto.push_back("A12B19C_hP64_194_ab2fk_efh2k_d;3;64;194;5;13;hP64;a,c/a,z4,z5,z6,z7,x8,x9,z9,x10,z10,x11,z11;-;Fe12O19Pb;Magnetoplumbite (PbFe12O19)");
    vproto.push_back("ABC4D_cP28_198_a_a_ab_a;4;28;198;6;8;cP28;a,x1,x2,x3,x4,x5,y5,z5;S6_{5};AlNaO4Si;alpha-Carnegieite (NaAlSiO4, S65)");
    vproto.push_back("AB2C4D_cP32_198_a_2a_ab_a;4;32;198;6;9;cP32;a,x1,x2,x3,x4,x5,x6,y6,z6;S6_{6};CaNa2O4Si;Na2CaSiO4 (S66)");
    vproto.push_back("A2B4C_cP56_198_ab_2a2b_2a;3;56;198;5;15;cP56;a,x1,x2,x3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8;-;Cu2O4Se;Cubic Cu2OSeO3");
    vproto.push_back("AB2_cI36_199_b_c;2;18;199;4;5;cI36;a,x1,x2,y2,z2;C26_{a};NO2;C26a (NO2) (obsolete)");
    vproto.push_back("A3B11C3_cP68_201_be_efh_g;3;68;201;5;8;cP68;a,x2,x3,x4,x5,x6,y6,z6;-;Bi3O11Ru3;Bi3Ru3O11");
    vproto.push_back("AB3C6D12_cF88_202_a_bc_e_h;4;22;202;6;4;cF88;a,x4,y5,z5;J2_{4};CoK3N6O12;K3Co(NO2)6 (J24)");
    vproto.push_back("A4BC12_cI34_204_c_a_g;3;17;204;5;3;cI34;a,y3,z3;-;Fe4LaP12;LaFe4P12");
    vproto.push_back("AB2_cI36_204_d_g;2;18;204;4;4;cI36;a,x1,y2,z2;C26;NO2;NO2 (Modern, C26)");
    vproto.push_back("A7BC12_cI40_204_bc_a_g;3;20;204;5;3;cI40;a,y4,z4;-;Mn7NaO12;NaMn7O12");
    vproto.push_back("A6BC_cP32_205_d_b_a;3;32;205;5;4;cP32;a,x3,y3,z3;-;F6NaSb;NaSbF6");
    vproto.push_back("A2B6C_cP36_205_c_d_a;3;36;205;5;5;cP36;a,x2,x3,y3,z3;G2_{1};N2O6Pb;Pb(NO3)2 (G21)");
    vproto.push_back("A4B_cP40_205_cd_c;2;40;205;4;6;cP40;a,x1,x2,x3,y3,z3;D1_{1};I4Sn;SnI4 (D11)");
    vproto.push_back("A7B2C_cP40_205_bd_c_a;3;40;205;5;5;cP40;a,x3,x4,y4,z4;K6_{1};P2O7Zr;ZrP2O7 High-Temperature (K61)");
    vproto.push_back("A2B6C6D_cP60_205_c_d_d_a;4;60;205;6;8;cP60;a,x2,x3,y3,z3,x4,y4,z4;J1_{10};Br2(H2O)6O6Zn;Zn(BrO3)2.6H2O (J110)");
    vproto.push_back("A2B6CD6_cP60_205_c_d_a_d;4;60;205;6;8;cP60;a,x2,x3,y3,z3,x4,y4,z4;H6_{4};N2(NH3)6NiO6;H64 [Ni(NO3)2(NH3)6] (obsolete)");
    vproto.push_back("A2BC4_cP84_205_d_ac_2d;3;84;205;5;11;cP84;a,x2,x3,y3,z3,x4,y4,z4,x5,y5,z5;-;B2CaO4;CaB2O4 (IV)");
    vproto.push_back("AB12CD8E2_cP96_205_a_2d_b_cd_c;5;96;205;7;12;cP96;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7;-;Cr(H2O)12NaO8S2;NaCr(SO4)2.12H2O Alum");
    vproto.push_back("AB24CD20E2_cP192_205_a_4d_b_c3d_c;5;192;205;7;24;cP192;a,x3,x4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11;H4_{15};AlH24NaO20S2;gamma-Alum [AlNa(SO4)2.12H2O, H415]");
    vproto.push_back("AB24CD28E2_cP224_205_a_4d_b_2c4d_c;5;224;205;7;28;cP224;a,x3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13;H4_{13};AlH24KO20S2;alpha-Alum [KAl(SO4)2.12H2O, H413]");
    vproto.push_back("AB2C36D2E20F2_cP252_205_a_c_6d_c_c3d_c;6;252;205;8;32;cP252;a,x2,x3,x4,x5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13,x14,y14,z14;H4_{14};AlC2H36N2O20S2;beta-Alum [Al(NH3CH3)2(SO4)2.12H2O, H414]");
    vproto.push_back("A2B3_cP60_212_bcd_ace;2;60;212;4;7;cP60;a,x3,x4,y5,x6,y6,z6;D5_{7};Fe2O3;Maghemite (gamma-Fe2O3, D57)");
    vproto.push_back("A3B2_cP20_213_d_c;2;20;213;4;3;cP20;a,x1,y2;-;Mg3Ru2;Mg3Ru2");
    vproto.push_back("A2BC3_cP24_213_c_a_d;3;24;213;5;3;cP24;a,x2,y3;-;Al2CMo3;Al2Mo3C");
    vproto.push_back("ABC4_cF24_216_b_a_e;3;6;216;5;2;cF24;a,x3;H0_{5};ClKO4;High-Temperature Cubic KClO4 (H05)");
    vproto.push_back("AB_cF40_216_de_ce;2;10;216;4;3;cF40;a,x3,x4;-;AlN;AlN (cF40)");
    vproto.push_back("AB4C8_cF52_216_a_e_2e;3;13;216;5;4;cF52;a,x2,x3,x4;-;GaMo4S8;GaMo4S8");
    vproto.push_back("A13BC18D20E5_cF228_216_dh_b_fh_2eh_ce;5;57;216;7;11;cF228;a,x4,x5,x6,x7,x8,z8,x9,z9,x10,z10;S0_{8};Al13ClF18O20Si5;Zunyite [Al13(OH,F)18Si5O20Cl, S08]");
    vproto.push_back("A16B40C12D6E5_cF316_216_eh_e2g2h_h_f_be;5;79;216;7;15;cF316;a,x2,x3,x4,x5,x6,x7,x8,z8,x9,z9,x10,z10,x11,z11;-;F16O40Ti12Y6Zn5;Murataite [(Y,Na)6(Zn,Fe)5Ti12O29(O,F)10F4]");
    vproto.push_back("A45B11_cF448_216_bd4efg5h_ac2eh;2;112;216;4;21;cF448;a,x5,x6,x7,x8,x9,x10,x11,x12,x13,z13,x14,z14,x15,z15,x16,z16,x17,z17,x18,z18;-;Cd45Sm11;Sm11Cd45");
    vproto.push_back("AB_cI16_217_c_c;2;8;217;4;3;cI16;a,x1,x2;-;AlN;AlN (cI16)");
    vproto.push_back("A4B24C13_cI82_217_c_deg_ag;3;41;217;5;7;cI82;a,x2,x4,x5,z5,x6,z6;-;As4Cu12S13;Tennantite (Cu12As4S13)");
    vproto.push_back("A3BC4D12E3_cP46_218_d_a_e_i_c;5;46;218;7;5;cP46;a,x4,x5,y5,z5;S6_{2};Al3ClNa4O12Si3;Sodalite [Na4(AlSiO4)3Cl, S62]");
    vproto.push_back("A3B4C4D4E16F4G3_cP76_218_c_e_e_e_ei_e_d;7;76;218;9;9;cP76;a,x3,x4,x5,x6,x7,x8,y8,z8;S6_{9};Al6Ca2.4K1.6Na4O30S1.5Si6;Hauyne [(Na0.5Ca0.3K0.2)8(Al6Si6O24)(SO4)1.5, S69]");
    vproto.push_back("AB_cI24_220_a_b;2;12;220;4;1;cI24;a;-;AlN;AlN (cI24)");
    vproto.push_back("A4B12C3_cI76_220_c_e_a;3;38;220;5;5;cI76;a,x2,x3,y3,z3;S1_{5};Bi4O12Si3;Eulytine (Bi4(SiO4)3, S15)");
    vproto.push_back("A7B12C19_cI152_220_bc_2d_ace;3;76;220;5;8;cI152;a,x3,x4,x5,x6,x7,y7,z7;K7_{4};Al14Ca12O33;Mayenite (12CaO.7Al2O3, K74, C12A7)");
    vproto.push_back("AB9C3_cI208_220_c_3e_e;3;104;220;5;14;cI208;a,x1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5;G5_{2};AlO9P3;Al(PO3)3 (G52)");
    vproto.push_back("AB12C_cP14_221_a_h_b;3;14;221;5;2;cP14;a,x3;-;CaH12Y;Predicted High-Pressure YCaH12");
    vproto.push_back("A3B40CD12_cP112_224_d_e3k_a_k;4;112;224;6;10;cP112;a,x3,x4,z4,x5,z5,x6,z6,x7,z7;-;(H3O)3O40PW12;H3PW12O40.3H2O");
    vproto.push_back("A5B40CD12_cP116_224_cd_e3k_a_k;4;116;224;6;10;cP116;a,x4,x5,z5,x6,z6,x7,z7,x8,z8;H4_{16};(H2.6O)5O40PW12;12-phosphotungstic acid [H3PW12O40.5H2O (H416)]");
    vproto.push_back("A27B52CD12_cP184_224_dl_eh3k_a_k;4;184;224;6;14;cP184;a,x3,x4,x5,z5,x6,z6,x7,z7,x8,z8,x9,y9,z9;-;H15O46PW12;Dodecatungstophosphoric Acid Hexahydrate [H3PW12O40.6H2O]");
    vproto.push_back("A2BC6D_cF40_225_c_a_e_b;4;10;225;6;2;cF40;a,x4;-;Ba2MnO6W;Double Perovskite (Ba2MnWO6)");
    vproto.push_back("A10B_cF44_225_cf_b;2;11;225;4;2;cF44;a,x3;-;H10La;LaH10 High-Tc Superconductor");
    vproto.push_back("A9B8_cF68_225_af_ce;2;17;225;4;3;cF68;a,x3,x4;D8_{9};Co9S8;Co9S8 (D89)");
    vproto.push_back("ABC6D8E2_cF72_225_b_a_e_f_c;5;18;225;7;3;cF72;a,x4,x5;H5_{8};ClFNa6O8S2;Sulphohalite [Na6ClF(SO4)2, H58]");
    vproto.push_back("A6B9CD2E6_cF96_225_e_bf_a_c_e;5;24;225;7;4;cF96;a,x4,x5,x6;J2_{5};C12Cu3Fe2(H2O)xN12;Cu3[Fe(CN)6]2.xH2O (J25, x approx 3)");
    vproto.push_back("AB30C16D3_cF200_225_a_ej_2f_bc;4;50;225;6;6;cF200;a,x4,x5,x6,y7,z7;J2_{1};AlF6H12N3;(NH4)3AlF6 (J21)");
    vproto.push_back("AB_cF32_227_c_d;2;8;227;4;1;cF32;a;L1_{3} (I);CuPt;Cubic CuPt (L13 (I), D4)");
    vproto.push_back("A2B_cF96_227_abf_cd;2;24;227;4;2;cF96;a,x5;D6_{2};O2Sb;D62 (Sb2O4) (obsolete)");
    vproto.push_back("A11B4_cF120_227_acdf_e;2;30;227;4;3;cF120;a,x4,x5;-;Ga2O3;gamma-Ga2O3");
    vproto.push_back("A16B2C_cF152_227_eg_d_a;3;38;227;5;4;cF152;a,x3,x4,z4;-;H16Li2Mg;Predicted Li2MgH16 High-Temperature Superconductor (250GPa)");
    vproto.push_back("A18B2C3_cF184_227_fg_d_ac;3;46;227;5;4;cF184;a,x4,x5,z5;-;Al18Cr2Mg3;Mg3Cr2Al18");
    vproto.push_back("A22B_cF184_227_cdfg_a;2;46;227;4;4;cF184;a,x4,x5,z5;-;Zn22Zr;Zn22Zr");
    vproto.push_back("A2BCD3E6_cF208_227_e_c_d_f_g;5;52;227;7;5;cF208;a,x3,x4,x5,z5;G7_{3};C2ClMgNa3O6;G73 [Northupite, Na3MgCl(CO3)2] (obsolete)");
    vproto.push_back("A4B2C6D16E_cF232_227_e_d_f_eg_a;5;58;227;7;6;cF232;a,x3,x4,x5,x6,z6;H5_{6};C4Mg2Na6O16S;H56 [Tychite, Na6Mg2SO4(CO3)4)] (obsolete)");
    vproto.push_back("A29B40CD12_cF656_227_ae2fg_e3g_b_g;4;164;227;6;15;cF656;a,x3,x4,x5,x6,x7,z7,x8,z8,x9,z9,x10,z10,x11,z11;H4_{21};(H2O)29O40PW12;H3PW12O40.29H2O (H421)");
    vproto.push_back("A21B_cI44_229_bdh_a;2;22;229;4;2;cI44;a,y4;B23;AgI;alpha-AgI (B23)");
    vproto.push_back("A2B3C12D12_cI232_230_a_c_h_h;4;116;230;6;7;cI232;a,x3,y3,z3,x4,y4,z4;J2_{3};Al2Ca3H12O12;Ca3Al2(OH)12 (J23)");
    // -------------------------------------------------------------------------
    // carbo-nitride prototypes (from DX) //DX20210728
    // -------------------------------------------------------------------------
    // quaternaries
    vproto.push_back("A2BCD4_mP16_11_2e_e_e_4e;4;16;11;6;20;mP16;a,b/a,c/a,beta,x1,z1,x2,z2,x3,z3,x4,z4,x5,z5,x6,z6,x7,z7,x8,z8;-;C2La1Li1N4;C2La1Li1N4 (ICSD #420440)");
    vproto.push_back("AB2C2D_mC24_15_a_f_f_e;4;12;15;6;11;mC24;a,b/a,c/a,beta,y2,x3,y3,z3,x4,y4,z4;-;Ag1C2N2Na1;Ag1C2N2Na1 (ICSD #65697)");
    vproto.push_back("A6BC6D2_hP15_147_g_a_g_d;4;15;147;6;9;hP15;a,c/a,z2,x3,y3,z3,x4,y4,z4;-;C6Fe1N6Pb2;C6Fe1N6Pb2 (ICSD #51493)");
    vproto.push_back("A6B2C6D_hP15_147_g_d_g_a;4;15;147;6;9;hP15;a,c/a,z2,x3,y3,z3,x4,y4,z4;-;C6Cd2N6Os1;C6Cd2N6Os1 (ICSD #417821)");
    vproto.push_back("A6B2CD6_hP15_147_g_d_a_g;4;15;147;6;9;hP15;a,c/a,z2,x3,y3,z3,x4,y4,z4;-;C6Cd2Fe1N6;C6Cd2Fe1N6 (ICSD #417826)");
    vproto.push_back("A3B6CD6_hP16_162_g_k_a_k;4;16;162;6;6;hP16;a,c/a,x3,z3,x4,z4;-;Ag3C6Co1N6;Ag3C6Co1N6 (ICSD #16959)");
    vproto.push_back("A6B2CD6_cF60_216_f_bc_a_f;4;15;216;6;3;cF60;a,x4,x5;-;C6Cu2Fe1N6;C6Cu2Fe1N6 (ICSD #28653)");
    vproto.push_back("A6BC6D_cF56_225_e_a_e_b;4;14;225;6;3;cF56;a,x3,x4;-;C6Cd1N6Pd1;C6Cd1N6Pd1 (ICSD #6093)");
    vproto.push_back("A3BCD3_cF64_225_e_c_ab_e;4;16;225;6;3;cF64;a,x4,x5;-;C6Cs2Mn2N6;C6Cs2Mn2N6 (ICSD #248046)");
#endif
    // done now produce

    // FROM PROTO LIST
    for(uint i=0;i<vproto.size();i++) {
      vproto_label.push_back("");
      vproto_nspecies.push_back(0);
      vproto_natoms.push_back(0);
      vproto_spacegroup.push_back(0);
      vproto_nunderscores.push_back(0);
      vproto_nparameters.push_back(0);
      vproto_Pearson_symbol.push_back("");
      vproto_params.push_back("");
      vproto_Strukturbericht.push_back("");
      vproto_prototype.push_back("");
      vproto_dialect.push_back("");

      anrl::vproto2tokens(vproto.at(i),
          vproto_label.at(i),vproto_nspecies.at(i),vproto_natoms.at(i),vproto_spacegroup.at(i),
          vproto_nunderscores.at(i),vproto_nparameters.at(i),vproto_Pearson_symbol.at(i),
          vproto_params.at(i),vproto_Strukturbericht.at(i),vproto_prototype.at(i),vproto_dialect.at(i));
    }

    return vproto.size();
  }
}

// ***************************************************************************
namespace anrl { 
  vector<string> getANRLParameters(string _anrl_label, string library, int choice, bool keep_original_lattice_parameter){
    //choice indicates choice of parameter set
    //all: grabs all
    //"": looks for one
    //part1: set from part1
    //part2: set from part2
    //part1:0
    //part1:1
    //part2:0
    //part2:1
    //...

    stringstream message;

    vector<string> tokens;
    vector<string> vparameters;

    string anrl_label = _anrl_label;

    // check if number suffix, i.e., predefined structure, with atomic volume scaling
    string number_id = "";
    if(aurostd::substring2bool(anrl_label,"-")){
      aurostd::string2tokens(_anrl_label,tokens,"-");
      anrl_label = tokens[0];
      number_id = tokens[1];
    }

    // add option to keep original scaling
    //bool keep_original_lattice_parameter = false;

    // A WORD TO THE WISE: When adding new preset parameter values for a given 
    // label, add them to the end.  Otherwise, you will mess up the connection
    // to the previously defined preset labels (e.g., -001, -002, etc.)

    // ---------------------------------------------------------------------------
    // Part 1
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="part1"){
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_aP12_1_4a_8a"){
        vparameters.push_back("5.417,1.0,1.0,90.0,90.0,90.0,0.001,0.002,0.003,0.4966,0.0001,0.5036,0.5001,0.502,0.0011,-0.0006,0.5013,0.5038,0.3857,0.3832,0.384,0.1149,0.6114,0.8846,0.8854,0.1157,0.6143,0.6153,0.8865,0.1141,0.6151,0.6132,0.6137,0.8854,0.3818,0.1149,0.1147,0.8856,0.3841,0.3857,0.1161,0.8842");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_aP16_1_4a_4a_8a"){
        vparameters.push_back("6.554,1.00061031431,1.92662496186,100.43475,100.46074,107.53,0.3267,0.582,0.177,0.565,-0.0132,0.4424,0.5217,0.3883,0.6767,-0.0744,0.6254,-0.0574,0.0338,0.0476,0.2599,0.0831,0.6072,0.4974,-0.0131,0.0949,0.7583,0.5449,0.1443,-0.0022,-0.0211,0.5213,0.2073,0.2907,0.5956,-0.0183,-0.0616,0.0602,0.4998,0.5068,-0.0175,0.2448,0.4596,0.0397,0.708,0.5326,0.352,0.4818,0.0,0.0,0.0,-0.078,0.569,0.7448");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_aP6_2_2i_i"){
        vparameters.push_back("4.56,1.54824561404,1.62280701754,80.2,106.96667,98.2,0.557,0.73,0.165,0.82,0.803,0.695,0.397,0.639,0.463");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_aP4_2_aci"){
        vparameters.push_back("3.307,2.24130631993,0.844572119746,89.06,85.15,85.7,0.572,0.259,0.433");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_3_bc3e_2e"){
        vparameters.push_back("4.1605,0.992524936907,1.78370388174,101.3752,0.15907,0.73859,0.02399,0.752,0.18927,0.38562,0.71473,0.64074,0.48963,0.20196,0.18802,0.18244,0.0,0.69651,0.38098,0.58564,0.17797");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP4_4_2a"){
        vparameters.push_back("3.104,2.42042525773,1.53350515464,92.71,0.25,0.23,0.48,0.48,0.0,0.02");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC12_5_3c"){
        vparameters.push_back("7.42,0.578167115903,1.90026954178,92.0,0.05,0.27,0.245,0.63,0.3,0.4,0.245,0.43,0.07");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_mC10_8_ab_a_a"){
        vparameters.push_back("5.72204,0.9978207073,0.722908263486,90.498,0.5515,-0.0994,0.0,0.0,0.523,0.4492,0.288,0.2434,0.3729");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC144_9_24a_12a"){
        vparameters.push_back("18.524,0.270092852516,1.28535953358,105.82,0.5749,0.351,0.8182,0.0707,0.34,0.8476,0.7315,0.138,0.4851,0.2509,0.144,0.5152,0.4155,0.352,0.6741,-0.0873,0.352,0.6434,0.8773,0.164,-0.0787,0.416,0.168,-0.0639,0.7741,0.145,0.7538,0.2336,0.143,0.7402,0.6195,0.341,0.5847,0.0811,0.343,0.5661,-0.0034,0.011,0.6062,0.3533,0.489,0.5665,0.6498,0.005,0.6711,0.1524,0.496,0.7805,0.8636,0.499,0.7328,0.3361,0.003,0.8333,0.0052,0.493,0.7398,0.1369,0.011,-0.0732,0.4927,0.492,0.8868,0.5,0.468,0.5,0.2252,0.491,0.5898,0.2744,0.021,-0.0845,0.0507,0.041,0.5642,0.2036,0.447,0.7347,-0.0802,0.049,0.6225,0.5751,0.043,0.7955,0.4247,0.048,0.6971,0.2643,0.444,0.5386,0.8023,0.449,0.7661,0.6453,0.041,0.6027,0.8531,0.463,-0.0984,0.4493,0.466,-0.0642,0.2244,0.059,-0.0395,0.0697,0.049,0.8702");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP4_11_e_e"){
        vparameters.push_back("2.8837,1.42393452856,1.61854561848,82.062,0.0387,0.8252,0.5887,0.7184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mP10_11_e_e_ef"){
        vparameters.push_back("4.63,1.20259179266,1.52203023758,110.21,0.121,0.1745,0.3531,0.7086,0.4009,0.1165,0.8544,0.5361,0.6943");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP16_11_8e"){
        vparameters.push_back("6.183,0.779880316998,1.77308749798,101.79,0.345,0.162,0.767,0.168,0.128,0.34,0.657,0.457,0.025,0.618,0.473,0.653,0.328,-0.074,0.869,0.894");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC6_12_a_i"){
        vparameters.push_back("7.189,0.613019891501,0.705105021561,90.04,0.6879,0.2889");
        vparameters.push_back("4.47,1.2212528,0.78724832,71.55,0.476,0.72");  // 002, binary metal-oxide prototype (ICSD #48214)
        vparameters.push_back("4.8754,0.57720392,0.98500636,70.4457,0.262,0.231");  // 003, binary metal-oxide prototype (ICSD #88720)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC34_12_ah3i2j"){
        vparameters.push_back("11.93871,0.876392843113,0.658278825769,129.00411,0.22,0.854,0.241,0.663,0.745,0.566,0.238,0.355,0.232,-0.037,0.333,0.35,0.586");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_12_g_ij"){
        vparameters.push_back("5.914,1.73047007102,1.03956712885,108.25,0.1662,0.2147,0.2263,0.2518,0.32131,0.2248");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_mC14_12_a2i_i"){
        vparameters.push_back("9.188,0.430343926861,0.705158902917,97.56,0.14286,0.42857,0.28571,0.85714,0.42857,0.28571");
        vparameters.push_back("5.2154,0.57732868,1.4994056,116.8579,0.333,0.5,0.628,0.942,0.85,0.275");  // 002, binary metal-boride prototype (ICSD #20326)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC4_12_i"){
        vparameters.push_back("5.403,0.635387747548,0.940033314825,132.32,0.106,0.173");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mP12_13_e_a_2g"){
        vparameters.push_back("8.95,0.500335195531,1.63360893855,145.35,0.5182,0.2986,0.0278,0.0003,0.2821,0.4045,0.2366");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP84_13_21g"){
        vparameters.push_back("9.21,0.99348534202,2.45385450597,106.1,0.30089,0.20127,0.18147,0.17387,0.03262,0.11695,0.05014,-0.05231,0.18035,-0.07589,0.78099,0.11634,0.79463,0.67872,0.1738,0.68463,0.51532,0.10402,0.56601,0.44932,0.17224,0.42424,0.27741,0.11672,0.0412,0.39067,0.07245,-0.00092,0.15881,0.04497,0.78847,0.13878,0.07346,0.7486,-0.09081,0.04464,0.53574,0.87264,0.06842,0.50833,0.63715,0.03304,0.30515,0.63715,0.06617,0.25041,0.40555,0.0442,0.146,0.38905,0.17219,0.86038,0.10055,0.17357,0.59606,0.82384,0.1694,0.41856,0.64581,0.16732,-0.05418,0.32296,0.2006");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_14_2e_e"){
        vparameters.push_back("5.1505,1.01186292593,1.03238520532,99.23,0.07,0.3317,0.3447,0.4496,0.7569,0.4792,0.2754,0.0395,0.2083");
        vparameters.push_back("5.4972,0.880103,1,120.47,0.8881,0.2186,0.7666,0.61,0.7024,0.7014,0.7722,0.9898,0.9889"); // Friedrich-binary-oxide (ICSD #8217) 
        vparameters.push_back("5.7038,0.78679126,0.93574459,122.6,0.106,0.212,0.209,0.401,0.703,0.299,0.239,0.979,0.026"); // Friedrich-binary-oxide (ICSD #647610) 
        vparameters.push_back("5.4128,0.90027343,0.90766332,123.7682,0.192,0.404,0.885,0.321,0.594,0.145,0.235,0,0.01");  // 004, binary metal-nitride prototype (ICSD #240759)
        vparameters.push_back("6.8392,0.7690958,0.78517955,131.2879,0.95,0.22,0.41,0.538,0.75,0.998,0.25,0,0.5");  // 005, binary metal-oxide prototype (ICSD #15983)
        vparameters.push_back("6.0818,0.79729685,0.77871683,131.75,0.056,0.653,0.774,0.425,0.273,0.962,0.309,0.955,0.091");  // 006, binary metal-oxide prototype (ICSD #154035)
        vparameters.push_back("6.023,0.79146605,0.78333057,132.0431,0.109,0.618,0.848,0.435,0.232,0.947,0.276,0.963,0.064");  // 007, binary metal-oxide prototype (ICSD #154036)
        vparameters.push_back("5.6891,0.83573852,0.97018861,121.453,0.106,0.189,0.193,0.392,0.712,0.268,0.265,0.1005,0.985");  // 008, binary metal-oxide prototype (ICSD #173153)
        vparameters.push_back("6.6815,0.76779166,0.79173838,131.3679,0.071,0.667,0.726,0.45,0.241,0.971,0.275,0.959,0.967");  // 009, binary metal-oxide prototype (ICSD #659226)
        vparameters.push_back("13.2594,0.31675641,0.84845468,141.1937,0.085,0.662,0.93,0.165,0.837,0.07,0.375,0.324,0.25");  // 010, binary metal-carbide prototype (ICSD #24074)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP32_14_8e"){
        vparameters.push_back("9.31,0.866809881847,1.38023630505,93.13333,0.437,0.185,0.084,0.246,0.273,-0.023,0.24,0.102,0.828,0.05,-0.08,0.852,0.157,0.669,-0.09,0.142,0.66,0.09,0.368,0.746,0.16,0.334,0.021,0.21");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP64_14_16e"){
        vparameters.push_back("15.018,0.979691037422,0.585231056066,93.61,0.18313,0.14063,0.03451,0.22856,0.28408,0.12262,0.35548,0.31907,-0.00548,0.47826,0.28776,0.16131,0.52853,0.14438,0.09345,0.47966,0.04033,0.27102,0.36296,-0.02818,0.15123,0.22521,0.04261,0.2343,0.09552,0.48601,0.14213,0.01298,0.58883,0.27815,-0.01931,0.71476,0.12135,0.08347,0.82945,0.18553,0.19177,0.81338,0.00963,0.3102,0.73961,0.14402,0.30834,0.59137,0.04778,0.24353,0.50553,0.23353");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5_mC28_15_f_e2f"){
        vparameters.push_back("12.786,0.387533239481,0.427968090099,97.03333,0.5727,0.106,0.311,0.077,0.0958,0.0952,0.4213,0.7127,0.0726,0.3138");
        vparameters.push_back("13.332,0.378173,0.988741,155.097,0.3306,0.892,0.5091,0.2798,0.8281,0.6701,0.3752,0.4039,0.3353,0.051"); // Friedrich-binary-oxide (ICSD #1422) , BA decoration
        vparameters.push_back("12.74,0.38328100471,0.436491365777,105.02,0.6465,0.3598,0.25874,0.2512,0.1088,0.4421,0.4709,0.2055,0.0723,0.1256");  // 003, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_15_c_e"){
        vparameters.push_back("4.6837,0.730747058949,1.0950317057,120.34,0.4184");
        vparameters.push_back("5.85,0.59487179,0.94017094,72.5,0.5");  // 002, binary metal-oxide prototype (ICSD #27659)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC48_15_ae3f_2f"){
        vparameters.push_back("7.1356,1.73344918437,1.00532541062,120.34,0.1163,0.266,0.1234,0.9401,0.3114,0.1038,0.3282,0.0172,0.2117,0.4782,0.14033,0.10833,0.07227,0.50682,0.15799,0.54077");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6D2_mC40_15_e_e_3f_f"){
        vparameters.push_back("9.79,0.901123595506,0.548518896834,105.81,0.3082,-0.0942,0.3888,0.4123,0.8659,0.1365,0.2411,0.6799,0.1468,0.4802,0.0124,0.2117,0.4057,0.7764"); //DX20210428 - equivalent to diopside (CaMg(SiO3)2, S4_{1}, http://aflow.org/prototype-encyclopedia/ABC6D2_mC40_15_e_e_3f_f.S4_1.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oP12_16_ag_cd_2u"){
        vparameters.push_back("5.61,1.01069518717,1.61319073084,0.2,0.26,0.125,0.74,0.8,0.63");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_18_ab_3c"){
        vparameters.push_back("8.32,1.15865384615,0.579326923077,0.0,0.0,0.25,0.25,0.0,0.25,0.5,0.5,0.124,0.309,0.382");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_19_2a_a"){
        vparameters.push_back("7.764,0.909582689335,0.558088614116,0.185,0.07,0.465,0.055,0.765,-0.008,0.884,-0.011,0.391");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC24_20_abc_c"){
        vparameters.push_back("8.74,0.576659038902,0.942791762014,0.3336,0.4403,0.2453,0.1971,0.2713,0.33154,0.03589,0.81143");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP2_25_b_a"){
        vparameters.push_back("2.8102,1.87104120703,1.0769696107,0.0,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP24_28_acd_2c3d"){
        vparameters.push_back("16.54,0.533252720677,0.269649334946,0.0,0.319,0.014,0.018,0.042,0.617,0.042,0.624,0.334,0.5,0.503,0.301,0.042,0.632,0.636,0.5,0.619,0.036,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP16_31_a_ab_2ab"){
        vparameters.push_back("7.43,0.869448183042,0.831763122476,0.8268,0.0,0.1514,0.4983,0.8226,0.6454,0.1436,0.1166,0.2466,0.3255,-0.0134,0.2598,0.3364,0.6184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_33_a_a"){
        vparameters.push_back("5.2857,1.11007056776,0.659950432298,0.1996,0.5867,0.2506,0.002,0.2003,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP32_33_a_3a_4a"){
        vparameters.push_back("9.11,1.01866081229,1.1613611416,0.2187,0.4807,0.2031,0.4418,0.2052,0.0015,0.4488,0.1967,0.4146,0.1422,0.9176,0.2246,0.191,0.2506,0.2228,0.3424,0.5361,0.0415,0.0069,0.5876,0.2212,0.3355,0.546,0.3761");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_36_2a_a"){
        vparameters.push_back("4.624,1.46820934256,2.69139273356,0.333,0.0,0.061,0.134,0.395,0.366");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC8_38_e_a_b"){
        vparameters.push_back("3.875,1.17470967742,1.59019354839,0.0,0.6144,0.155,0.2914");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_38_de_ab"){
        //DX20201105 - changed y4 to keep in SG#38 (see comments in part 1) vparameters.push_back("4.684,1.81084543126,1.0269000854,0.06,0.5,0.17,0.56,0.17,0.0");
        vparameters.push_back("4.684,1.81084543126,1.0269000854,0.06,0.5,0.17,0.56,0.1701,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_oC20_41_a_2b"){
        vparameters.push_back("6.388,1.00485284909,1.7778647464,0.0,0.673,0.327,0.376,0.827,0.673,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC24_41_2a_2b"){
        vparameters.push_back("6.478,1.0,1.87635072553,0.01,0.238,0.342,0.158,0.125,0.25,0.25,-0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oF72_43_ab_3b"){
        vparameters.push_back("11.66,1.91595197256,0.58833619211,0.0,0.125,0.13889,0.0,0.02222,0.08056,0.18333,0.15278,-0.01389,-0.18333,0.0625,0.125,0.27778");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI4_44_a_b"){
        vparameters.push_back("4.92,0.973577235772,0.535569105691,0.0,0.425");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C7D_oP13_47_t_aq_eqrs_h"){
        vparameters.push_back("3.8187,1.01691675177,3.05567339671,0.3554,0.1579,0.3771,0.3788,0.18445");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_51_e_f"){
        vparameters.push_back("4.7549,0.661969757513,1.0209678437,0.8125,0.3125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP20_56_ce_e"){
        vparameters.push_back("4.911,2.53797597231,1.10201588271,0.029,0.147,0.058,0.861,0.044,0.128,0.179");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_oP16_57_d_c_d_d"){
        vparameters.push_back("6.707,0.997614432682,1.13627553303,0.208,0.7704,0.2871,0.889,0.4154,0.605,0.1087");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_57_d_d"){
        vparameters.push_back("6.09556,0.900425883758,0.850291031505,0.8593,0.0628,0.255,0.0096");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP6_58_a_g"){
        vparameters.push_back("6.24,1.03044871795,0.673076923077,0.275,0.325");
        vparameters.push_back("4.704,0.917942176871,0.601615646259,0.66667,0.25");
        vparameters.push_back("4.4446,1.22049228277,0.761913333033,0.2004,0.3787");
        vparameters.push_back("4.025,1.1965466,0.67115528,0.8682,0.5945");  // 004, binary metal-nitride prototype (ICSD #160620)
        vparameters.push_back("5.013,1.2288051,0.58567724,0.414,0.068");  // 005, binary metal-nitride prototype (ICSD #166465)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_59_a_b"){
        vparameters.push_back("3.15,1.29841269841,2.20634920635,0.051,0.277");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP6_59_a_a_a"){
        vparameters.push_back("5.68,0.700704225352,1.01056338028,0.1499,0.4237,0.6255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP8_59_bf_a"){
        vparameters.push_back("5.162,0.842115459124,0.877760557923,0.67125,0.329,0.505,0.174");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP16_61_c_c"){
        vparameters.push_back("6.471,1.27538247566,1.31757070005,0.136,0.072,0.108,0.456,0.119,0.872");
        vparameters.push_back("6.19,1.0193861,1.1437803,0.91,0.0327,0.9517,0.7614,0.0857,0.8729");  // 002, binary carbo-nitride prototype (ICSD #15870)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_61_2c_c"){
        vparameters.push_back("9.174,0.375953782429,0.560061042075,0.0095,0.1491,0.1835,0.2314,0.111,0.5366,0.1289,0.0972,0.8628");
        vparameters.push_back("15.035,0.363418689724,0.372929830396,0.028,0.634,0.171,0.168,0.221,0.081,0.1181,0.0252,0.3378");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP20_62_3c_2c"){
        vparameters.push_back("11.282,0.339443361106,0.994947704308,0.2922,0.19181,0.4504,0.877,0.6246,0.5611,-0.02937,0.17398,0.64939,-0.03603");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oP20_62_c_cd_a"){
        vparameters.push_back("5.4224,1.41099881971,0.996661994689,0.4877,-0.0084,0.0313,0.0586,0.288,0.537,0.213");
        vparameters.push_back("6.2152,1.4341775,1.0315356,0.9878,0.4505,0.3956,0.5401,0.6976,0.0554,0.2"); //Friedrich-ternary-oxide (ICSD #97463)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oP20_62_2cd_c"){
        vparameters.push_back("5.464,0.810395314788,1.36749633968,0.22451,0.65626,0.55801,0.6466,0.05131,0.36362,0.13079,0.0579,0.06543");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_62_c_2c_c"){
        vparameters.push_back("6.018,0.630741110003,2.4086075108,0.2522,0.8276,0.6221,0.095,0.8706,0.8244,0.226,0.06333");
        vparameters.push_back("12.41,0.31933924,0.43424658,0.3765,0.884,0.3273,0.0819,0.4236,0.6848,0.1309,0.8856");  // 002, ternary metal-carbo-nitride prototype (ICSD #75040)
        vparameters.push_back("5.5566,0.69605514,2.111903,0.8857,0.5978,0.667,0.5591,0.0745,0.6391,0.60706,0.3644");  // 003, ternary metal-carbo-nitride prototype (ICSD #410915)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_62_2c_c"){
        vparameters.push_back("4.918,0.7600650671,1.44550630338,0.038,0.282,0.674,0.562,0.202,0.611");
        vparameters.push_back("12.735,0.468237141735,0.339615233608,0.733,0.125,0.508,0.722,0.874,0.447");
        vparameters.push_back("7.6204,0.595008136056,1.1869718125,0.125,0.4217,0.0202,0.837,0.2377,0.0959");
        vparameters.push_back("3.875,1.64232258065,1.89496774194,0.004,0.758,0.24,0.07,0.24,0.39"); //DX20201019 - moved from Part 2 to here
        vparameters.push_back("9.579,0.29710826,0.43334377,0.174,0.545,0.094,0.259,0.891,0.1");  // 005, binary metal-nitride prototype (ICSD #187446)
        vparameters.push_back("5.52,0.62862319,1.1780797,0.65,0.55,0.04,0.33,0.738,0.891");  // 006, binary metal-oxide prototype (ICSD #56696)
        vparameters.push_back("5.22,0.64750958,1.3250958,0.628,0.6035,0.021,0.357,0.739,0.917");  // 007, binary metal-oxide prototype (ICSD #181282)
        vparameters.push_back("5.181,0.58733835,1.1713955,0.574,0.641,0.341,0.024,0.884,0.756");  // 008, binary metal-oxide prototype (ICSD #182577)
        vparameters.push_back("6.3708,0.60259308,1.1028285,0.8538,0.9201,0.976,0.3349,0.72666,0.60119");  // 009, binary metal-oxide prototype (ICSD #380398)
        vparameters.push_back("11.42,0.376532399299,0.805604203152,0.103,0.119,0.614,0.842,0.811,0.108");  // 010, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_62_c_c"){
        vparameters.push_back("10.42,0.349328214971,0.411708253359,0.375,0.333,0.139,0.389");
        vparameters.push_back("5.24160,0.606723137973,1.12622100122,0.0056,0.1952,0.1879,0.5696");
        vparameters.push_back("5.495,0.536123748863,0.737579617834,0.125,0.69,-0.18,0.125");
        vparameters.push_back("11.18,0.356171735242,0.387209302326,0.3507,0.0201,0.61937,0.3806");
        vparameters.push_back("5.454,0.609644297763,1.10542720939,0.2005,0.5741,0.0058,0.1993"); //DX20210428 - equivalent to eta-NiSi (B_{d}, http://aflow.org/prototype-encyclopedia/AB_oP8_62_c_c.NiSi.html, part 3)
        vparameters.push_back("6.9613,0.83485843,0.5325873,0.1136,0.7544,0.3592,0.4045");
        vparameters.push_back("6.12,0.5,0.74509804,0.029,0.397,0.177,0.877");  // 007, binary metal-boride prototype (ICSD #24701)
        vparameters.push_back("6.4,0.46875,0.875,0.5,0.48,0.26,0.76");  // 008, binary metal-boride prototype (ICSD #153291)
        vparameters.push_back("5.251,0.57912779,0.75261855,0.037,0.36,0.18,0.875");  // 009, binary metal-boride prototype (ICSD #612863)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_62_c_cd"){
        vparameters.push_back("5.09,1.3257367387,0.888605108055,0.39,0.05,0.036,0.852,0.186,0.063,0.328");
        vparameters.push_back("7.09,1.0603667,0.67616361,0.376,0.943,0.026,0.085,0.838,0.928,0.384");  // 002, binary metal-nitride prototype (ICSD #165990)
        vparameters.push_back("5.08,1.4409449,0.89566929,0.716,0.64,0.038,0.198,0.763,0.944,0.344");  // 003, binary metal-carbide prototype (ICSD #167667)
        vparameters.push_back("5.287,1.238623,0.89249102,0.36928773,0.55653649,0.03201707,0.33818756,0.82208818,0.9404134,0.15557256");  // 004, binary metal-carbide prototype (ICSD #617486)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_oP40_62_cd_3c2d"){
        vparameters.push_back("4.526,1.54882898807,2.68272205038,0.4594,0.5629,0.0579,0.6261,0.2501,0.2063,0.2619,0.4165,0.0288,0.0291,0.3428,0.0565,0.0642,0.8119,0.2509,0.0657,0.0218");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oP8_62_2c"){
        vparameters.push_back("6.663,0.708839861924,0.73345339937,0.464,0.292,0.181,0.658");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_63_c_2c_c"){
        vparameters.push_back("3.577,4.56863293263,1.09538719597,0.06109,-0.0558,0.1792,0.33096");
        vparameters.push_back("3.5814,5.0438934,1.3735132,0.95974,0.60233,0.78501,0.19337");  // 002, ternary metal-nitride prototype (ICSD #85528)
        vparameters.push_back("4.4247,3.7013583,0.98573915,0.697,0.5638,0.2498,0.9044");  // 003, metallic ternary prototype (ICSD #240096)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_63_2c_c"){
        vparameters.push_back("3.73,3.94638069705,0.983914209115,0.061,0.75,0.396");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_63_c_c"){
        vparameters.push_back("2.9782,2.64253575985,0.985360284736,0.436,0.14525");
        vparameters.push_back("1.0.0185797325,1.41421356236,0.707106781182,0.625,0.125"); // Lederer-13
        vparameters.push_back("3.31,1.4441088,1.1389728,0.3174,0.6888");  // 003, binary metal-nitride prototype (ICSD #163951)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oC4_63_c"){
        vparameters.push_back("2.8444,2.06331739558,1.73379271551,0.10228");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oC8_64_f"){
        vparameters.push_back("4.523,1.69378730931,1.0002210922,0.1549,0.081");
        vparameters.push_back("3.3136,3.16211974891,1.32070859488,0.10168,0.08056");
        vparameters.push_back("7.11906,0.654575182679,1.37596817557,0.15485,0.1175");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC80_64_efg_efg_df"){
        vparameters.push_back("10.922,0.866233290606,0.682933528658,0.84657,0.0946,0.9271,0.5886,0.276,-0.0792,0.2314,0.27981,-0.0113,0.1278,0.3415,0.2438,0.1245,0.175,0.2231");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_65_j_g"){
        vparameters.push_back("5.971,1.1314687657,0.468263272484,0.28,0.22");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oC16_65_ah_bej"){
        vparameters.push_back("8.031,0.926410160628,0.491595069107,0.25,0.225");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC8_65_a_bf"){
        vparameters.push_back("5.82068,1.35259626023,0.493507631411");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oF8_69_a_b"){
        vparameters.push_back("6.08,0.903782894737,0.851973684211");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oF8_70_a"){
        vparameters.push_back("3.1587,1.82613100326,3.21714629436");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oF24_70_e_a"){
        vparameters.push_back("8.2671,0.580614725841,1.0342804611,0.4615");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oF128_70_4h"){
        vparameters.push_back("10.4646,1.22947843205,2.33988876785,0.14415,0.04732,0.0486,0.29277,0.2269,0.25406,0.21598,0.28022,0.32618,0.21405,0.15761,0.37947");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oI6_71_a_i"){
        vparameters.push_back("3.144,0.994910941476,2.44179389313,0.339");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oI6_71_a_g"){
        vparameters.push_back("2.75984,2.9999963766,1.4241115427,0.35333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oI12_72_j_a"){
        vparameters.push_back("9.583,0.585829072316,0.578837524783,0.1182,0.2088");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tI12_82_c_g_a"){
        vparameters.push_back("4.3404,1.53216293429,0.256,0.2566,0.3722");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI14_82_bc_a_g"){
        vparameters.push_back("5.55,1.85585585586,0.26,0.25,0.13");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP16_84_cej_k"){
        vparameters.push_back("6.429,1.02830922383,0.46779,0.25713,0.19361,0.30754,0.22904");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5_tI18_87_h_ah"){
        vparameters.push_back("10.164,0.37111373475,0.2797,-0.0589,0.3752,0.6856");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_tI10_87_a_h"){
        vparameters.push_back("5.72,0.623076923077,0.4,0.8");
        vparameters.push_back("5.188,0.41094834,0.722,0.423");  // 002, binary metal-carbide prototype (ICSD #187144)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP12_92_b_a"){
        vparameters.push_back("4.957,1.39001412144,0.3047,0.2381,0.1109,0.1826");
        vparameters.push_back("4.7985,1.584391,0.036776,0.153876,0.270876,0.195876"); // Friedrich-binary-oxide (ICSD #161691) //DX20210428 - equivalent to alpha-O2Te (paratellurite, http://aflow.org/prototype-encyclopedia/A2B_tP12_92_b_a.TeO2.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP36_96_3b_ab"){
        vparameters.push_back("7.464,1.15487674169,0.41,0.445,0.132,0.4,0.117,0.123,0.296,0.344,0.297,0.143,0.326,0.12,0.248");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP12_96_ab"){
        vparameters.push_back("5.51889,1.25999974633,0.0849,0.1752,0.3792,0.2742");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tP5_99_bc_a_b"){
        vparameters.push_back("4.046,1.02308452793,0.0,0.8973,0.4517,0.3785");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP8_113_a_ce"){
        vparameters.push_back("6.871,0.606622034638,0.206,0.1797,0.476");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4D_tI16_121_d_a_i_b"){
        vparameters.push_back("5.46,1.96428571429,0.245,0.132");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_122_a_b_d"){
        vparameters.push_back("5.289,1.97069389299,0.2574");
        vparameters.push_back("4.11,1.8399027,0.82");  // 002, ternary metal-carbo-nitride prototype (ICSD #44110)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C_tP7_123_b_ci_a"){
        vparameters.push_back("4.207,1.61516520086,0.312");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP4_123_a_ce"){
        vparameters.push_back("4.158,0.864357864358");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_a_d"){
        vparameters.push_back("2.8,1.31071428571");
        vparameters.push_back("4.44,0.643243243243");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_d_a_f"){
        vparameters.push_back("3.8611,0.828649866618");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP10_127_g_ah"){
        vparameters.push_back("7.3364,0.530232811733,0.3841,0.1821");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_tP8_129_c_b_a_c"){
        vparameters.push_back("3.6736,2.60540069686,0.6793,0.2246");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP4_129_ac"){
        vparameters.push_back("4.897,0.69185215438,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_c_a_c"){
        vparameters.push_back("4.11,1.76301703163,0.6497,0.2058");
        vparameters.push_back("4.002,1.5527236,0.3638,0.7938");  // 002, ternary metal-oxide prototype (ICSD #32743)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_129_ac_c"){
        vparameters.push_back("4.0006,1.52584612308,0.27,0.7");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_129_a_c"){
        vparameters.push_back("3.9645,1.26008323874,0.2368");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_129_c_c"){
        vparameters.push_back("3.107,1.90505310589,0.1,0.65");
        vparameters.push_back("1.0,2.82842712472,0.875,0.375"); // Lederer-42
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_131_c_e"){
        vparameters.push_back("4.9073,1.24500234345");
        vparameters.push_back("2.979,1.7582746"); //Friedrich-binary-oxide (ICSD #26598)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP50_134_b2m2n"){
        vparameters.push_back("8.74,0.575514874142,0.0048,0.1685,0.1305,0.628,0.1695,0.5228,0.1635,0.0753,0.3383,0.1485");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="sigma_tP30_136_bf2ij" || anrl_label == "A_tP30_136_bf2ij"){ //DX20210630 - combined disordered and unary version (they are structurally the same)
        vparameters.push_back("10.59,0.532011331445,0.1033,0.3667,0.0383,0.5608,0.2354,0.3183,0.27"); //DX20210630 - the explicit sigma parameters are slightly different: 8.7966,0.518177,0.10136,0.36878,0.03651,0.56609,0.23933,0.31733,0.25202
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP8_136_g_f"){
        vparameters.push_back("4.75,0.576842105263,0.31,0.336");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_136_f_a"){
        vparameters.push_back("4.5922,0.644005052045,0.30496");
        vparameters.push_back("4.86,0.56995885,0.305");  // 002, binary metal-oxide prototype (ICSD #647647)
        vparameters.push_back("3.9342,1.2762442,0.6092");  // 003, binary metal-carbide prototype (ICSD #88057)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP4_136_f"){
        vparameters.push_back("3.957,1.29112964367,0.098");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP16_138_j"){
        vparameters.push_back("8.56,0.714953271028,0.375,-0.083,0.857");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI16_139_cde_e"){
        vparameters.push_back("3.9993,4.3215062636,0.37498,0.11886");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI4_139_e"){
        vparameters.push_back("3.34916,1.94217355994,0.819");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_tI14_139_a_e_ce"){
        vparameters.push_back("3.7817,3.50337149959,0.36075,0.1824");
        vparameters.push_back("4.6113,3.24412,0.645,0.848"); //Friedrich-ternary-oxide (ICSD #157402), equivalent to K2NiF4 (http://aflow.org/prototype-encyclopedia/A4B2C_tI14_139_ce_e_a.html, part 3)
        vparameters.push_back("3.629,3.4169193,0.387,0.16");  // 003, ternary metal-oxide prototype (ICSD #86753)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_tI26_139_fij_a"){
        vparameters.push_back("8.47,0.584415584416,0.361,0.278");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI2_139_a"){
        //DX 20191218 [this is the a' parameter (fct) vs the a parameter (bct)] vparameters.push_back("4.6002,1.07523585931");
        vparameters.push_back("3.25283,1.52061313499"); //DX 20191218 [CORRECT PARAMETERS]
        vparameters.push_back("3.932,0.823499491353");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI8_139_h"){
        vparameters.push_back("4.33184,0.574102459925,0.17916");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI8_139_bd_a"){
        vparameters.push_back("3.8537,2.22744375535");
        vparameters.push_back("1.0,2.0"); // Lederer-46
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI6_139_a_e"){
        vparameters.push_back("3.2064,2.44754241517,0.3353");
        vparameters.push_back("1.0,4.24264068707,0.3333333333"); // Lederer-44
        vparameters.push_back("4.4371,1.7910121,0.596"); // Friedrich-binary-oxide (ICSD #24248)
        vparameters.push_back("4.7583,1.66103,0.4025"); // Friedrich-binary-oxide (ICSD #38245) //DX20210428 - equivalent to C2Ca-I (C11_{a}, http://aflow.org/prototype-encyclopedia/A2B_tI6_139_e_a.html, part 3)
        vparameters.push_back("5.01,1.1816367,0.61");  // 005, binary metal-oxide prototype (ICSD #20275)
        vparameters.push_back("5.384,1.2706166,0.6089");  // 006, binary metal-oxide prototype (ICSD #24729)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5_tI18_139_i_ah"){
        vparameters.push_back("8.91,0.361391694725,0.328,0.348");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tI10_139_de_a"){
        vparameters.push_back("4.53,2.45033112583,0.38");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B_tI18_139_hi_a"){
        vparameters.push_back("8.312,0.468840230991,0.333,0.327");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI6_139_d_a"){
        vparameters.push_back("4.1,1.22682926829");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_140_h_a"){
        vparameters.push_back("6.04,0.804635761589,0.158");
        vparameters.push_back("4.4988,1.1500622,0.6605");  // 002, binary metal-carbide prototype (ICSD #181489)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_b_ah"){
        vparameters.push_back("6.017,1.44241316271,0.231");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI16_140_ab_h"){
        vparameters.push_back("8.03,0.87297633873,0.179");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_tI24_141_h_b_a"){
        vparameters.push_back("6.6042,0.905423821205,0.066,0.1951");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI4_141_a"){
        vparameters.push_back("5.8318,0.545611989437");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_tI28_141_ad_h"){
        vparameters.push_back("5.765,1.63781439722,0.0278,0.2589");
        vparameters.push_back("8.16,1.1568627,0.485,0.735");  // 002, binary metal-oxide prototype (ICSD #643198)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_141_e_a"){
        vparameters.push_back("3.785,2.51360634082,0.08306");
        vparameters.push_back("4.126,3.47697527872,0.2915"); //DX20201103 - moved from part 2
        vparameters.push_back("4.2,1.9047619,0.617");  // 003, binary metal-nitride prototype (ICSD #30593)
        vparameters.push_back("3.78,2.515873,0.6975");  // 004, binary metal-oxide prototype (ICSD #200392)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI16_141_e_e"){
        vparameters.push_back("3.108,5.45045045045,0.227,0.071");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI24_141_2e_e"){
        vparameters.push_back("4.046,6.28917449333,0.125,0.289,-0.051");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI8_141_a_b"){
        vparameters.push_back("3.325,3.42255639098");
        vparameters.push_back("1.0,1.99999999997"); // Lederer-53
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tI80_141_ceh_3h"){
        vparameters.push_back("7.5937,4.26037373086,0.2044,0.5201,0.3324,0.516,0.2547,0.494,0.0859,0.4667,0.4164");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tI96_142_e_ab_2g"){
        vparameters.push_back("10.914,1.77396005131,0.0375,0.2482,0.3197,-0.0867,0.0923,0.1117,0.0025");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_147_g_ad"){
        vparameters.push_back("7.636,0.369264012572,0.25,0.33333,0.0,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR16_148_cf_cf"){
        vparameters.push_back("6.29713,1.8633345667,0.11546,0.21,0.10706,0.81289,0.19519,0.1848,0.6754,0.3468");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR8_148_c_f"){
        vparameters.push_back("7.49626,2.75900649124,0.33333,0.088,0.755,0.421");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR26_148_b2f_a2f"){
        vparameters.push_back("15.659,0.335334312536,0.054,0.346,0.098,0.754,0.15699,0.6,0.555,0.84401,0.599,0.252,0.65501,0.098");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR10_148_c_f_c"){
        vparameters.push_back("5.0884,2.76815894977,0.35537,0.1464,0.22174,0.56249,0.95095");
        vparameters.push_back("5.2482146,2.8404053,0.342,0.156,0.97,0.26,0.54");  // 002, ternary metal-oxide prototype (ICSD #15989)
        vparameters.push_back("5.183002,2.6774052,0.1439,0.3507,0.2029,0.5575,0.989");  // 003, ternary metal-oxide prototype (ICSD #153674)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_150_ef_bd"){
        vparameters.push_back("5.85,0.589743589744,0.875,0.26,0.6");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_151_3c_2a"){
        vparameters.push_back("6.017,2.87518697025,0.8889,0.5556,0.8889,0.1111,0.0731,0.5556,0.4444,0.0731,0.2222,0.77778,0.0731");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_152_c_a"){
        vparameters.push_back("4.914,1.10012210012,0.4699,0.413,0.2668,0.214");
        vparameters.push_back("5.291,1.1591382,0.453,0.408,0.303,0.21503333");  // 002, binary metal-oxide prototype (ICSD #41493)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP3_152_a"){
        vparameters.push_back("4.3662,1.13453346159,0.2254");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_154_a_b"){
        vparameters.push_back("4.145,2.29095295537,0.7198,0.4889");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR8_155_c_de"){
        vparameters.push_back("4.91608,2.53341483458,0.237,0.43,0.07");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hR5_155_e_c"){
        vparameters.push_back("5.73296,1.24097324942,0.2521,0.2449");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR6_160_b_b"){
        vparameters.push_back("9.619,0.327466472606,0.00019,0.26362,0.7288,0.39161");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR6_160_3a_3a"){
        vparameters.push_back("3.01791,7.34847294982,0.0,0.22222,0.77778,0.08333,0.30556,0.86111");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR10_161_a_a_b"){
        vparameters.push_back("5.2542,2.64091583876,0.2875,0.0128,0.74643,0.14093,0.36263");
        vparameters.push_back("5.485751,2.5028469,0,0.7419,0.4088,0.5408,0.9728");  // 002, ternary metal-oxide prototype (ICSD #9645)
        vparameters.push_back("5.148249,2.6927581,0.775,0,0.6395,0.8964,0.306");  // 003, ternary metal-oxide prototype (ICSD #28296)
        vparameters.push_back("5.5357687,2.4202479,0,0.7492,0.8091,0.2648,0.6945");  // 004, ternary metal-oxide prototype (ICSD #51036)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_162_ad_k"){
        vparameters.push_back("4.917,0.929021761237,0.325,0.272");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2CD2_hP36_163_h_i_bf_i"){
        vparameters.push_back("7.384,2.37716684724,0.01,0.833,0.33333,0.03833,0.141,0.03167,0.365,0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP5_164_ad_d"){
        vparameters.push_back("4.0282,1.21409066084,0.648,0.149");
        vparameters.push_back("4.1971,1.2861976,0.31689,0.76948");  // 002, binary metal-nitride prototype (ICSD #162795)
        vparameters.push_back("3.8572,1.7764181,0.33541,0.76913");  // 003, binary metal-nitride prototype (ICSD #162796)
        vparameters.push_back("3.938,1.2496191,0.606,0.231");  // 004, binary metal-nitride prototype (ICSD #169727)
        vparameters.push_back("3.3888,1.670916,0.3521,0.752");  // 005, binary metal-oxide prototype (ICSD #160203)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_164_a_d"){
        vparameters.push_back("4.24,1.61320754717,0.252");
        vparameters.push_back("3.387,1.9790375,0.82221");  // 002, binary metal-carbide prototype (ICSD #280743)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_165_adg_f"){
        vparameters.push_back("6.308,1.03994927077,0.167,0.666,0.356,0.028,0.096");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR2_166_a_b"){
        vparameters.push_back("3.13,4.78594249201");
        vparameters.push_back("3.0357007,2.486313");  // 002, binary metal-oxide prototype (ICSD #82236)
        vparameters.push_back("9.88622,1.00556532224");  // 003, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR2_166_c"){
        vparameters.push_back("3.7595,2.7815666977,0.22754");
        vparameters.push_back("2.456,4.08957654723,0.16667");
        vparameters.push_back("3.289,3.42991790818,0.0543");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR1_166_a"){
        vparameters.push_back("5.07846,0.968139947937");
        vparameters.push_back("3.45741,1.92728082582");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B6_hR13_166_ah_3c"){
        vparameters.push_back("4.757,5.4319949548,0.167,0.346,0.448,0.09,0.59001");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR3_166_ac"){
        vparameters.push_back("3.62036,7.25049442597,0.22222");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hR5_166_c_ac"){
        vparameters.push_back("4.36914,6.96313919902,0.399,0.208");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hR7_166_a2c_c"){
        vparameters.push_back("3.011,6.9511790103,0.186,0.33333,0.075");
        vparameters.push_back("3.0113138,6.9505254,0.167,0.314,0.425");  // 002, binary metal-boride prototype (ICSD #20326)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR12_166_2h"){
        vparameters.push_back("4.908,2.56022616137,0.0104,0.65729,0.2206,0.6323");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR4_166_a_b_c"){
        vparameters.push_back("3.5561,5.44557239673,0.2667");
        vparameters.push_back("1.0,4.89897948553,0.25"); //Lederer-59
        vparameters.push_back("3.2903652,5.0765915,0.2375"); //Friedrich-ternary-oxide (ICSD #28288)
        vparameters.push_back("2.8496891,5.6601982,0.899"); //Friedrich-ternary-oxide (ICSD #246912) //DX20210428 - equivalent to rhombohedral delafossite (CuFeO2, http://aflow.org/prototype-encyclopedia/ABC2_hR4_166_a_b_c.CuFeO2.html, part 3)
        vparameters.push_back("2.9850154,6.2010049,0.8881");  // 005, ternary metal-oxide prototype (ICSD #4149)
        vparameters.push_back("3.8300036,4.4647567,0.108");  // 006, ternary metal-oxide prototype (ICSD #18102)
        vparameters.push_back("3.6300006,4.7024827,0.892");  // 007, ternary metal-oxide prototype (ICSD #18106)
        vparameters.push_back("3.3042225,5.2820657,0.8945");  // 008, ternary metal-oxide prototype (ICSD #55687)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hR105_166_bc9h4i"){
        vparameters.push_back("10.96,2.17974452555,0.3848,0.3843,0.21309,0.4895,0.21780,0.3873,0.56899,0.1991,0.50609,0.1983,0.68740,0.1032,0.49209,0.9933,0.66980,0.1008,0.83740,0.0025,0.16801,0.3622,0.58109,0.0976,0.3765,0.68261,0.2024,0.1673,0.55209,0.8921,0.1777,0.3473,0.0033");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_hR7_166_g_a"){
        vparameters.push_back("4.33304,3.13251204697,0.16667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR10_167_a_b_e"){
        vparameters.push_back("5.285,2.62039735099,0.85756666666667");
        vparameters.push_back("4.988,3.42040898156,0.5067");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hR10_167_c_e"){
        vparameters.push_back("4.7607,2.72957758313,0.35216,0.5561");
        vparameters.push_back("4.95878,2.82762,0.654,0.9381"); // Friedrich-binary-oxide (ICSD #94768)
        vparameters.push_back("4.7600055,2.7289929,0.315,0.553");  // 003, binary metal-oxide prototype (ICSD #43732)
        vparameters.push_back("4.7948913,2.7319666,0.64726,0.4167");  // 004, binary metal-oxide prototype (ICSD #99783)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP18_180_fi_bd"){
        vparameters.push_back("5.198,2.54136206233,0.163,0.1141");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_180_d_j"){
        vparameters.push_back("4.42758,1.43826876081,0.16559");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_180_j_c"){
        vparameters.push_back("4.9977,1.09252256038,0.2072");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP8_182_c_g"){
        vparameters.push_back("4.8507,0.866967654153,0.3249");
        vparameters.push_back("4.7309,0.93863747,0.69255333");  // 002, binary metal-boride prototype (ICSD #184958)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_186_ab"){
        vparameters.push_back("2.47,2.75303643725,0.0,0.07143");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_186_ab_ab"){
        vparameters.push_back("3.08051,3.27374363336,0.18784,0.0,0.43671,0.24982");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_186_b_b"){
        vparameters.push_back("3.8227,1.63776911607,0.3748,0.0");
        vparameters.push_back("3.143,1.823099,0,0.375");  // 002, binary metal-nitride prototype (ICSD #153888)
        vparameters.push_back("4.084,1.4299706,0.5,0.91");  // 003, binary metal-nitride prototype (ICSD #162195)
        vparameters.push_back("2.965,2.0856661,0,0.3748");  // 004, binary metal-nitride prototype (ICSD #168369)
        vparameters.push_back("3.925,1.16,0.375,0");  // 005, binary metal-carbide prototype (ICSD #189087)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_186_a2b_a2b"){
        vparameters.push_back("3.08129,4.90695780014,0.1254,0.0,0.29215,-0.0415,0.16675,0.8335");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C_hP18_186_2a3b_2ab_b"){
        vparameters.push_back("3.281,6.57726302956,0.155,0.345,0.0,0.248,0.045,0.261,0.455,0.367,0.137");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_186_b_a"){
        vparameters.push_back("2.51,2.66932270916,0.0,0.05");
        vparameters.push_back("3.172,1.6314628,0,0.38");  // 002, binary metal-nitride prototype (ICSD #159250)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP3_187_a_d_f"){
        vparameters.push_back("4.535,1.0769570011");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP2_187_d_a"){
        vparameters.push_back("2.9065,0.975950455875");
        vparameters.push_back("2.864,1.1075419");  // 002, binary metal-nitride prototype (ICSD #186245)
        vparameters.push_back("3.0026,0.81649237");  // 003, binary metal-oxide prototype (ICSD #166273)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_189_fg_bc"){
        vparameters.push_back("5.877,0.584822188191,0.256,0.589");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP6_191_a_h_b"){
        vparameters.push_back("3.04436,2.20489035462,0.2413");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5_hP6_191_a_cg"){
        vparameters.push_back("5.405,0.773913043478");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP1_191_a"){
        vparameters.push_back("3.2062,0.931195808122");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP4_191_bc_a"){
        vparameters.push_back("3.6576,1.05902777778");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_191_a_d"){
        vparameters.push_back("3.005,1.08276206323");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_191_h_e"){
        vparameters.push_back("4.237,1.71040830776,0.306,0.16");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_191_f_ad"){
        vparameters.push_back("5.279,0.806914188293");
        vparameters.push_back("5.1912,0.56069888");  // 002, binary metal-nitride prototype (ICSD #25659)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_194_ad_f"){
        vparameters.push_back("3.64,3.37362637363,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP6_194_h"){
        vparameters.push_back("4.40445,0.568892824303,0.44799");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_194_af_bf"){
        vparameters.push_back("3.01,4.85382059801,0.166,0.583");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_194_ac"){
        vparameters.push_back("3.77,3.2175066313");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP8_194_c_bf"){
        vparameters.push_back("5.088,1.76533018868,-0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_b_f"){
        vparameters.push_back("4.895,1.58324821246,0.045");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_194_c_d"){
        vparameters.push_back("2.50399,2.66023426611");
        vparameters.push_back("3.03,1.2508251");  // 002, binary metal-nitride prototype (ICSD #163950)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP8_194_d_a_f"){
        vparameters.push_back("2.86,4.48251748252,0.086");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_194_h_c"){
        vparameters.push_back("5.295,0.802077431539,0.8392");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_194_bc"){
        vparameters.push_back("2.464,2.72362012987");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_c_f"){
        vparameters.push_back("3.161,3.8895919013,0.6275");
        vparameters.push_back("2.7767,3.5081212,0.606");  // 002, binary metal-nitride prototype (ICSD #169882)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hP14_194_abdf_f"){
        vparameters.push_back("2.982,4.651240778,0.528,0.139");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP12_194_f_ah"){
        vparameters.push_back("5.223,1.64005360904,0.06286,0.83048");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_c_d_a"){
        vparameters.push_back("2.752,2.56468023256");
        vparameters.push_back("4.4075,1.3225184");  // 002, metallic ternary prototype (ICSD #57506)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP4_194_f"){
        vparameters.push_back("2.508,1.66786283892,0.05995");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_c_ad"){
        vparameters.push_back("4.186,1.22527472527");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_hP16_194_c_af_ef"){
        vparameters.push_back("2.988,7.82195448461,0.1543,0.605,0.0539");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP2_194_c"){
        vparameters.push_back("3.2093,1.62359393014");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP24_194_ef_fgh"){
        vparameters.push_back("4.824,3.28067993367,0.04598,0.84417,0.12514,0.16429");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_194_df_ce"){
        vparameters.push_back("3.976,4.12022132797,0.0637,0.10724");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_194_c_a"){
        vparameters.push_back("3.619,1.39375518099");
        vparameters.push_back("2.968,1.8696092");  // 002, binary metal-nitride prototype (ICSD #76384)
        vparameters.push_back("3.048,1.6135171");  // 003, binary metal-nitride prototype (ICSD #105123) //DX20210428 - equivalent to Fe2N (L'3_{0}, http://aflow.org/prototype-encyclopedia/AB_hP4_194_c_a.Fe2N.html, part 3)
        vparameters.push_back("2.7472,2.1177927");  // 004, binary metal-nitride prototype (ICSD #181300)
        vparameters.push_back("2.93,1.9931741");  // 005, binary metal-nitride prototype (ICSD #185559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP12_194_cg_f"){
        vparameters.push_back("5.052,1.63697545527,0.062");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cI40_197_cde_c"){
        vparameters.push_back("8.4295,0.1668,0.3345,0.6476,0.7484");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_cP12_198_a_a_a"){
        vparameters.push_back("5.881,-0.024,0.39,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP16_198_b_a"){
        vparameters.push_back("5.1305,0.2107,0.3689,0.2671,0.1159");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP8_198_2a"){
        vparameters.push_back("5.65,0.0699,-0.0378");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP8_198_a_a"){
        vparameters.push_back("5.63,-0.042,0.067");
        vparameters.push_back("4.48688,0.13652,0.8424");
        vparameters.push_back("4.46,0.558,0.822");  // 003, binary metal-carbide prototype (ICSD #168277)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cI16_199_a_a"){
        vparameters.push_back("6.3557,0.294,0.0347");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB32C48_cI162_204_a_2efg_2gh"){
        vparameters.push_back("14.16,0.8203,0.5998,0.1836,0.2942,0.8806,0.0908,0.8499,0.1748,0.6993,0.686,0.0969,0.332");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cI32_204_g_c"){
        vparameters.push_back("7.58,0.3431,0.8497");
        vparameters.push_back("7.4456,0.7702,0.7283");  // 002, binary metal-oxide prototype (ICSD #55462)
        vparameters.push_back("7.255,0.7087,0.2018");  // 003, binary metal-oxide prototype (ICSD #55469)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_cI26_204_g_a"){
        vparameters.push_back("7.58,0.184,0.691");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP8_205_c"){
        vparameters.push_back("5.65,0.05569");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP16_205_c_c"){
        vparameters.push_back("6.4162,0.1527,0.6297");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cP12_205_a_c"){
        vparameters.push_back("5.417,0.3851");
        vparameters.push_back("5.9431,0.427"); //DX20200703 - binary oxide for R. Friedrich (ICSD #87178)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C6_cI80_206_b_d_e"){
        vparameters.push_back("9.4,-0.0344,0.338,0.1,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI16_206_c"){
        vparameters.push_back("4.11971,0.1001");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP20_213_cd"){
        vparameters.push_back("6.315,0.06361,0.20224");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_cP8_215_d_e_a"){
        vparameters.push_back("5.3912,0.2372");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_cP5_215_a_e"){
        vparameters.push_back("3.878,0.265");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_cP8_215_a_c_e"){
        vparameters.push_back("5.28,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5_cF24_216_a_ce"){
        vparameters.push_back("6.1,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_cF12_216_b_c_a"){
        vparameters.push_back("6.24");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF8_216_c_a"){
        vparameters.push_back("5.4093");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cI10_217_c_a"){
        vparameters.push_back("5.45858,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI58_217_ac2g"){
        vparameters.push_back("8.911,0.31787,-0.08958,0.28194,0.64294,0.03457");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B8_cI52_217_ce_cg"){
        vparameters.push_back("8.8664,0.32774,0.10781,0.64421,0.68844,0.03674");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI16_220_c"){
        vparameters.push_back("5.2716,0.049");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_cI40_220_d_c"){
        vparameters.push_back("8.135,0.0492,0.2896");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP2_221_b_a"){
        vparameters.push_back("4.07925"); //DX20210428 - equivalent to NH4NO3 I (G0_{8}, http://aflow.org/prototype-encyclopedia/AB_cP2_221_a_b.NH4.NO3.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP6_221_c_d"){
        vparameters.push_back("4.2101");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_cP5_221_a_c_b"){
        vparameters.push_back("3.795");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB27CD3_cP32_221_a_dij_b_c"){
        vparameters.push_back("7.04,0.245,0.26");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_cP4_221_a_c"){
        vparameters.push_back("3.7402");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP1_221_a"){
        vparameters.push_back("3.34");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB11_cP36_221_c_agij"){
        vparameters.push_back("9.6,0.345,0.225,0.115");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB11CD3_cP16_221_a_dg_b_c"){
        vparameters.push_back("5.74,0.245");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP4_221_d_a"){
        vparameters.push_back("3.734");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_cP7_221_f_a"){
        vparameters.push_back("4.145,0.2117");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP8_223_c_a"){
        vparameters.push_back("4.556");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP46_223_dik"){
        vparameters.push_back("10.355,0.1837,0.1172,0.3077");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cP6_224_b_a"){
        vparameters.push_back("4.267");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B_cF32_225_bd_a"){
        vparameters.push_back("9.45"); //DX20210428 - equivalent to CuPt3 (partially occupied, L1_{a}, http://aflow.org/prototype-encyclopedia/AB7_cF32_225_b_ad.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_cF16_225_a_bc"){
        vparameters.push_back("5.853");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B16C7_cF128_225_acd_2f_be"){
        vparameters.push_back("11.48,0.25,0.875,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_cF52_225_i_a"){
        vparameters.push_back("7.468,0.666");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cF12_225_a_c"){
        vparameters.push_back("5.4631");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B23_cF116_225_e_acfh"){
        vparameters.push_back("10.65,0.2765,0.6191,0.6699");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cF16_225_a_c_b"){
        vparameters.push_back("5.95");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF4_225_a"){
        vparameters.push_back("3.61491");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB18C8_cF108_225_a_eh_f"){
        vparameters.push_back("10.56,0.325,0.65833,0.66");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF8_225_a_b"){
        vparameters.push_back("5.63931");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cF24_227_c_a"){
        vparameters.push_back("7.166");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cF96_227_e_cf"){
        vparameters.push_back("11.278,0.215,0.44");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF16_227_a_b"){
        vparameters.push_back("7.483");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF136_227_aeg"){
        vparameters.push_back("14.864,0.2624,0.1824,0.3701");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cF24_227_d_a"){
        vparameters.push_back("7.02");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF8_227_a"){
        vparameters.push_back("3.55");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_cF56_227_d_a_e"){
        vparameters.push_back("8.0832,0.7376");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cF48_227_c_e"){
        vparameters.push_back("8.6,0.245");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C3_cF112_227_c_de_f"){
        vparameters.push_back("11.087,0.7047,0.323");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cI2_229_a"){
        vparameters.push_back("3.155");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cI8_229_b_a"){
        vparameters.push_back("2.984");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cI14_229_c_b"){
        vparameters.push_back("6.226");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7_cI54_229_e_afh"){
        vparameters.push_back("11.618,0.6862,0.1704,0.6503");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB12C3_cI32_229_a_h_b"){
        vparameters.push_back("7.04,0.7625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C3_cI16_229_a_c_b"){
        vparameters.push_back("5.74");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cI112_230_af_g"){
        vparameters.push_back("11.411,0.0,0.625");
      }
    }
    // ---------------------------------------------------------------------------
    // Part 2
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="part2"){
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_aP6_2_aei_i"){
        vparameters.push_back("2.7804,1.00438785786,1.53100273342,78.28,76.53,70.42,0.259,0.415,0.663,0.192,0.183,0.255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B5_mP13_6_a7b_3a2b"){
        vparameters.push_back("6.5369054321,0.490874879533,1.43786810261,109.592,0.5119,0.7172,0.4546,0.4285,0.044,0.5911,0.37,0.0053,0.3887,0.2195,0.6277,0.0001,0.1818,0.4663,0.7456,0.5256,0.1826,0.7901,0.7824,0.2724,0.0,0.0,0.8188,0.8026,0.0123,0.2051");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP4_6_2b_2a"){
        vparameters.push_back("3.580975428,1.00027925162,1.00167550965,90.04,0.0,0.0,0.518,0.507,0.026,0.501,0.529,0.027");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_7_4a_2a"){
        vparameters.push_back("5.0942,0.627360527659,1.04603274312,90.38,0.498,-0.062,0.472,-0.023,0.574,0.143,0.777,-0.052,0.799,0.271,0.151,0.261,-0.001,0.808,0.36,0.494,0.649,0.658");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP18_7_6a_3a"){
        vparameters.push_back("6.5499839723,1.91328244277,1.22717557252,127.75,0.6753,0.4828,0.3393,0.4711,0.6498,0.3067,0.4174,0.3433,0.0893,0.2883,0.1551,0.3432,0.1472,0.1112,0.5582,0.0,0.0746,0.0,0.5693,0.07852,0.0933,0.1184,0.41567,0.2938,0.8473,0.27155,0.6656");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mP16_7_6a_2a"){
        vparameters.push_back("5.2778002048,0.97690325515,1.45210125431,91.762,0.5044,0.292,0.01,0.5764,0.215,0.586,0.0,0.209,0.0,0.0864,0.29,0.58,0.2874,0.0717,0.287,0.7924,0.4201,0.301,0.2874,0.014,0.0012,0.7994,0.528,0.078");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2_mP22_7_9a_2a"){
        vparameters.push_back("6.4165085102,0.999298672153,1.36910105356,93.39,0.4944,0.2369,0.0135,0.2965,0.2623,0.5568,0.4955,0.4453,0.2912,0.2817,0.0557,0.2805,0.8067,0.5491,0.0506,0.0,0.0309,0.0,0.6774,0.1539,0.7442,0.096,0.6343,0.3258,0.8638,0.2389,0.2526,0.6359,0.1217,0.4513,0.156,0.3761,0.11377");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_mC32_9_5a_3a"){
        vparameters.push_back("8.12077,0.718445418353,1.12797801194,115.809,0.009,-0.003,0.269,0.129,0.341,0.45,0.37,0.119,0.066,0.142,0.351,0.147,0.356,0.135,0.348,0.0,0.5182,0.0,0.136,0.2,0.309,0.365,0.2924,0.196");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_9_a_3a"){
        vparameters.push_back("3.5367,2.67777872028,0.986823875364,93.018,0.50691,0.35551,0.49997,0.28659,0.26502,0.27913,0.57076,0.06256,0.41996,0.42173,0.07037,0.56036");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_10_mn_bg"){
        vparameters.push_back("4.012,0.819541375872,2.93668993021,97.03,0.843,0.126,0.558,0.644");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mP16_10_mn_3m3n"){
        vparameters.push_back("3.678,0.71098423056,1.23817292007,91.0,0.263,0.339,0.08,0.059,0.795,0.341,0.438,-0.069,0.763,0.161,0.705,0.841,0.58,0.441,0.062,0.431");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mP8_10_ac_eh_mn"){
        vparameters.push_back("5.1177434753,0.86107854631,1.45123094959,90.021,0.6089,0.24179,0.1277,0.24913");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP6_10_en_am"){
        vparameters.push_back("5.1700416367,0.61508704062,1.49709864605,104.5,0.234,0.66,0.263,0.336");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mP8_10_2m2n"){
        vparameters.push_back("4.7302,0.527461840937,0.86332501797,106.1,0.1175,0.6746,0.5344,0.3333,0.1131,0.8977,0.4209,0.1319");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C2_mC22_12_aij_h_i"){
        vparameters.push_back("6.65,1.29563909774,0.704661654135,102.2,0.30503,0.38654,0.22171,0.22108,-0.08762,0.23655,0.15499,0.71826");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC16_12_4i"){
        vparameters.push_back("9.089,0.274617669711,0.451534822313,96.96,-0.0572,0.1206,0.4419,0.3467,0.7858,-0.0594,0.2715,0.4149");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_13_2g_ef"){
        vparameters.push_back("5.6255,0.61198115723,1.23780997245,127.44,0.1808,-0.004,0.155,0.346,0.225,0.345,0.273,0.573");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_14_e_a"){
        vparameters.push_back("5.5496,0.695689779444,1.15512829753,107.151,0.255,0.2573,0.3141");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B8_mP120_14_14e_16e"){
        vparameters.push_back("7.5889,0.766725085322,3.55545599494,106.136,0.25252,0.5751,0.3402,0.76924,0.888,0.47034,0.13489,0.3866,0.33222,0.8218,0.8146,0.42747,0.12506,0.2343,0.29191,0.78423,0.9464,0.38285,0.23276,0.2673,0.25885,0.69189,0.1533,0.38027,0.34844,0.4547,0.26591,0.63874,0.2281,0.42248,0.35819,0.6069,0.3061,0.67752,0.0969,0.46714,0.26777,0.7386,0.3843,0.8119,0.7461,0.51889,0.0602,0.3617,0.3547,0.8843,0.6723,0.4287,0.0438,0.1068,0.2871,0.822,0.8945,0.354,0.2273,0.1621,0.2315,0.6653,0.243,0.3497,0.4219,0.4796,0.2431,0.5754,0.3699,0.4209,0.4385,0.7352,0.3104,0.6409,0.1506,0.496,0.9409,0.7672,0.5381,0.7891,0.5835,0.5099,0.7334,0.7952,0.5403,0.2081,0.8842,0.3711,0.3975,0.7664,0.4019,0.2077,0.6717,0.4087");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_15_e_cf"){
        vparameters.push_back("3.282,2.64137720902,0.970749542962,91.9,0.86,0.577,0.068,0.169");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC24_15_2e2f"){
        vparameters.push_back("4.939,0.569143551326,0.838023891476,142.47,0.1012,0.3684,0.226,0.0672,0.2464,0.3443,0.1958,0.2227");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_17_abe_e"){
        vparameters.push_back("7.0499703347,1.11347517732,0.614184397158,0.893,0.878,0.379,0.225,0.522,0.202,0.275,0.022");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_19_a_3a"){
        vparameters.push_back("5.668,0.758997882851,0.528757939308,0.584,0.123,0.027,0.31,0.159,0.417,0.257,0.073,0.603,0.983,0.124,0.227");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC6_21_a_k"){
        vparameters.push_back("3.3982513,1.3943496174,1.40170688639,0.268");
        vparameters.push_back("3.343,1.73197726593,2.34519892312,0.34");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oF40_22_fi_ad_gh"){
        vparameters.push_back("6.4793068924,1.3977127159,1.54737394473,0.3134,0.3627,0.1144,0.0719");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oF8_22_a_c"){
        vparameters.push_back("5.5400291632,0.990433212996,0.937725631773");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oI32_23_ij2k_k"){
        vparameters.push_back("5.82463,1.24369101557,1.32254065924,0.04851,0.45153,0.75475,0.49255,0.20405,0.4421,0.23223,0.28616,0.76005,0.8216,0.36488");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2C12D2E_oI50_23_bcfk_i_3k_j_a"){
        vparameters.push_back("10.767,0.502554100492,1.4969815176,0.2511,0.3298,0.1693,0.2465,0.0107,0.1695,0.1308,0.2443,0.0826,0.3792,0.7558,0.0801,0.1294,0.7488,0.2546");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI16_23_ab_i_k"){
        vparameters.push_back("5.3999384419,1.15740740741,2.00555555555,0.28,0.25,0.2,0.115");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oI12_23_a_b_k"){
        vparameters.push_back("5.2501580231,1.06666666665,1.7219047619,0.21,0.2,0.115");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB7CD2_oI44_24_a_b3d_c_ac"){
        vparameters.push_back("7.0501914381,1.035035461,1.41546099291,0.7511,0.2496,0.1361,-0.0003,0.5002,0.7501,0.2213,0.3356,0.5662,0.0681,0.1362,-0.0648,0.0709,0.1385");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_26_abc_ab"){
        vparameters.push_back("4.6806,0.627034995513,1.05710806307,0.455,0.858,0.179,0.623,0.048,0.545,0.375,0.355,0.751,0.119,0.213");
        vparameters.push_back("5.0725824349,0.881353258932,1.48474034935,0.12,0.0,0.2484,0.461,0.246,0.289,0.3781,0.0862,0.253,0.652,0.12");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B_oP24_26_3a3b2c_ab"){
        vparameters.push_back("6.465013742,1.07113689096,1.87440061873,0.1628,0.2907,0.112,0.0,0.1452,0.8256,0.4859,0.1227,0.1146,0.0096,0.3286,0.1358,0.1565,0.2905,0.7214,0.2818,0.2489,0.2854,0.3947,0.249,0.0963,0.5436");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B4C16D_oP108_27_abcd4e_4e_16e_e"){
        vparameters.push_back("13.0409854622,0.999877309088,0.703110981603,0.0,0.005,0.03,0.05,0.373,0.379,0.296,0.122,0.125,0.249,0.369,0.129,0.28,0.122,0.376,0.26,0.058,0.751,0.449,0.237,0.574,0.086,0.257,0.012,0.554,0.483,0.247,0.028,0.729,0.183,0.412,0.162,0.684,0.26,0.23,0.164,0.679,0.651,0.299,0.6217,0.4,0.254,0.238,0.433,0.091,0.44,0.445,0.395,0.455,0.245,0.09,0.304,0.404,0.057,0.126,0.093,0.054,0.102,0.243,0.402,0.336,0.099,0.466,0.124,0.392,0.467,0.154,0.103,0.253,0.191,0.047,0.385,0.425,0.052,0.111,0.41,0.74634,0.2513,0.28218");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_29_2a_a"){
        vparameters.push_back("5.2594682584,0.963498098863,0.965209125484,0.639,0.068,0.0,0.771,0.537,0.106,0.53,0.267,0.356");
        vparameters.push_back("5.068,1.0378848,1.0017758,0.932,0.639,0.106,0.463,0.771,0,0.733,0.97,0.25");  // 002, binary metal-oxide prototype (ICSD #67004)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP12_29_a_2a"){
        vparameters.push_back("5.4179557747,1.0,1.0,0.5049,0.2419,0.0,0.615,0.135,0.3834,0.615,0.635,0.1134");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_29_a_a_a"){
        vparameters.push_back("5.2594682584,0.963498098863,0.965209125484,0.61885,0.63065,0.11668,0.50496,0.24091,0.0,0.61734,0.13129,0.37996");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C15_oP46_30_a2c_bc_a7c"){
        vparameters.push_back("5.4630061801,1.00183049606,3.84605528098,0.5,0.107,0.523,0.2442,0.513,0.019,0.3675,0.026,0.016,0.1074,0.003,0.02,0.1949,0.05,0.074,0.097,0.687,0.22,0.0844,0.254,0.226,0.413,0.52,0.105,0.505,0.718,0.31,0.305,0.763,0.271,0.695,0.752,0.289");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_30_2a_c_3c"){
        vparameters.push_back("4.480982758,1.71211783083,3.19772372237,0.8543,0.5,0.3011,0.721,0.4096,0.3607,0.7305,0.1733,0.5847,0.8629,0.0496,0.6302,0.8523,0.2991");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B2C2_oP34_32_a6c_c_c"){
        vparameters.push_back("8.4261784988,1.1450866366,0.701103726564,0.0,0.808,0.55,0.8946,0.614,0.7595,-0.0255,0.8197,0.8044,0.642,0.5367,0.6481,0.3894,0.7689,-0.0703,0.2917,0.8392,0.688,0.2793,0.67666,0.60773,0.0844,0.8643,0.8235,0.4147");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oP40_33_4a_6a"){
        vparameters.push_back("4.8437,1.71975968784,1.84873134174,0.6787,0.8416,0.0,0.1846,0.3432,0.7868,0.8115,0.6489,0.6972,0.6677,0.4696,0.9993,0.329,0.8313,0.8927,0.0248,0.4908,0.6292,0.4717,0.6647,0.6381,0.5145,0.6728,0.1212,0.8608,0.3301,0.8662,0.336,0.4992,0.9");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C_oP22_34_c_4c_a"){
        vparameters.push_back("6.2740008991,2.04032515143,1.38284985656,0.5,0.605,0.8135,0.499,0.7349,0.6796,0.009,0.743,-0.0925,0.3064,0.2388,0.5935,0.2196,0.7131,0.6459,0.513");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP6_34_a_c"){
        vparameters.push_back("5.8327827022,1.12083390481,0.548158688784,0.5,0.6881,0.8565,0.0097");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_oC22_35_a_ab3e_e"){
        vparameters.push_back("4.534758023,1.1408839779,2.72580110498,0.0,0.5838,-0.042,0.189,0.5461,0.0961,0.122,0.2982,0.1399,0.1866,0.0038");
        vparameters.push_back("3.62,5.359116,1.140884,0.9716,0.5554,0.9296,0.811,0.5177,0.9039,0.0936,0.2982,0.1115,0.8134,0.9754");  // 002, ternary metal-oxide prototype (ICSD #25378)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_36_a_a"){
        vparameters.push_back("5.825,0.945115879828,0.922403433476,0.25,0.0,0.081,0.83");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C2_oC36_37_d_c2d_d"){
        vparameters.push_back("5.8071602545,2.51110728429,0.821939039086,0.0,0.654,0.0584,0.0469,0.6705,0.0718,0.471,0.0932,0.1377,0.4004,0.1552,0.14836,0.0571");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oC40_39_2d_2c2d"){
        vparameters.push_back("7.47831,2.30292940517,0.749515599113,0.3449,0.5,0.0058,0.7654,0.173,0.5398,0.3889,0.5843,0.6154,0.3917,0.28104,0.16619,0.0206,0.10455,0.6073,0.0144");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9BC_oC44_39_3c3d_a_c"){
        vparameters.push_back("6.2123847492,1.90215711527,2.62509658728,0.5,0.2673,0.7493,0.6176,0.8879,0.6246,0.6059,0.6191,0.7476,0.2266,0.0857,0.2459,0.1789,0.5932,0.0649,0.18,0.5948,0.4281");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_40_a_2b_b"){
        vparameters.push_back("6.4580033647,1.33199132859,1.69634561784,0.0,0.3467,0.1773,0.8661,0.3301,0.7281,0.0093");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_40_b_3b"){
        vparameters.push_back("4.3850022361,5.92335058951,0.997331752147,0.83109,0.0,0.70417,0.0021,0.56971,0.4966,-0.07002,0.4978");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B3_oF52_42_2abce_ab"){
        vparameters.push_back("7.4494846573,1.70036689767,1.04688136975,0.76,0.27,0.0,0.31,0.06,0.79,0.6,0.17,0.11,0.07");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oF8_42_a_a"){
        vparameters.push_back("2.5000573158,1.33999999997,1.73599999999,0.0,0.333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oI20_45_c_b_c"){
        vparameters.push_back("5.9678340183,1.97721179624,0.981568364609,0.25,0.2729,0.629,0.493,0.2296,0.8629,0.474");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oI36_46_ac_bc_3b"){
        vparameters.push_back("6.9969353511,1.54780620265,0.898527940538,0.25,0.5253,0.0054,0.7207,0.7706,-0.0021,-0.0823,0.2996,0.7963,0.5295,0.6236,0.1199,0.006,0.3325,0.4952");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8CD_oP24_48_k_2m_d_b"){
        vparameters.push_back("6.3302356554,1.0,1.50710900473,0.4978,0.56,0.629,0.114,0.642,0.558,0.601");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_oP14_49_dehq_ab"){
        vparameters.push_back("3.6705354001,1.69539132804,2.12544314151,0.002,0.681");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C8D_oP24_49_g_q_2qr_e"){
        vparameters.push_back("6.3302356554,1.0,1.50710900473,0.5184,0.7996,0.258,-0.0697,0.3652,0.6307,0.2617,0.1943,0.8286");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oP28_50_ij_ac_ijm"){
        vparameters.push_back("5.5348069961,2.26684733515,0.987895212291,0.8726,0.445,0.3867,-0.076,0.5,0.743,0.226");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_oP48_50_3m_m_2m"){
        vparameters.push_back("11.0940438033,1.50045069408,0.472480620154,0.0515,0.8152,0.238,0.602,0.5567,0.315,0.604,0.8447,0.116,0.61286,0.66266,0.2344,0.11964,0.66839,0.2479,0.63183,0.00412,0.2439");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_52_2e_cd"){
        vparameters.push_back("7.2235008877,1.34576036547,1.32098013428,0.3175,0.6759,0.8271,0.1762,0.5576,0.5093,0.0419,0.8142");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP20_52_de_cd"){
        vparameters.push_back("15.5832022616,0.434585119311,0.426069989123,0.443,0.4294,0.0001,0.6539,0.064,0.0788");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_53_h_e_gh"){
        vparameters.push_back("11.6904618594,0.282132455361,0.78890717994,0.79514,0.3193,0.1198,0.3549,0.224,0.7493");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_53_e_g_hi"){
        vparameters.push_back("14.3630679002,0.312469539787,0.535821207264,0.6826,0.2856,0.3458,0.2708,0.6247,0.6057,0.3575");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_54_e_d_cf"){
        vparameters.push_back("5.3467489374,0.956627452436,1.81710213776,0.8667,0.6417,0.8902,0.4055,0.2686,0.5503");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_55_2g2h_gh"){
        vparameters.push_back("10.1600191136,1.45275590551,0.366929133855,0.0628,0.4022,0.1014,0.1118,0.2024,0.2667,0.226,0.0384,0.3532,0.2953,0.4192,0.1378");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oP16_55_ch_agh"){
        vparameters.push_back("5.4199981729,1.90405904061,0.730627306278,0.348,0.22,0.112,0.152,0.17,0.393");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oP16_55_2g2h"){
        vparameters.push_back("7.7886,0.613101199189,0.320442698303,0.6731,-0.037,0.8435,0.8087,-0.0454,0.8613,0.5704,0.8926");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP6_58_g_a"){
        vparameters.push_back("3.7572,2.89952624295,0.89063664431,0.3326,0.63309");
        vparameters.push_back("2.7,1.8259259,1.5296296,0.87249,0.40386");  // 002, binary metal-nitride prototype (ICSD #157283)
        vparameters.push_back("3.778,1.2916887,0.84912652,0.852,0.584");  // 003, binary metal-nitride prototype (ICSD #166463)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP6_59_a_b_a"){
        vparameters.push_back("3.301,1.14298697364,2.39612238716,0.32961,-0.04795,0.89243");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oP20_60_d_cd"){
        vparameters.push_back("8.4598035458,0.706855791961,0.724586288418,0.547,0.394,0.75,0.033,0.348,0.611,0.396");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP32_60_3d_d"){
        vparameters.push_back("7.3397836195,1.05524748967,1.03197678378,0.5016,0.7205,0.0322,0.2167,0.7591,0.2582,0.2197,0.5016,0.013,0.248,0.783,0.0291");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B8_oP120_60_7d_8d"){
        vparameters.push_back("13.85,0.824548736462,0.53357400722,0.382,0.3932,0.5783,0.2764,0.38,0.5411,0.2264,0.4631,0.4411,0.1289,0.4507,0.4069,0.0794,0.3554,0.4717,0.1276,0.2723,0.571,0.2251,0.2844,0.6054,0.2616,0.5329,0.393,0.094,0.5116,0.3344,0.0088,0.3466,0.4468,0.0917,0.2026,0.6185,0.2593,0.2232,0.6779,0.3972,0.4784,0.595,0.4192,0.3624,0.4724,0.4002,0.3492,0.6895");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP48_61_3c_3c"){
        vparameters.push_back("6.914,1.08128435059,1.38313566676,0.1297,0.5762,0.40803,0.1235,0.6328,0.54518,0.0057,0.4432,0.36289,0.2172,0.6275,0.346,0.2068,0.7225,0.5756,0.0095,0.4051,0.2704");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oP20_62_2c_3c"){
        vparameters.push_back("5.5329,0.511305102207,2.07339731425,0.1008,0.2055,0.2432,-0.0464,0.0157,0.4015,0.1808,0.7737,0.8691,-0.0688");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oP28_62_ac_2cd_c"){
        vparameters.push_back("10.193,0.586382811734,0.466202295693,0.2774,-0.0085,0.0913,0.7657,0.4474,0.2215,0.094,0.4262,0.1628,0.0331,0.2777");
        vparameters.push_back("9.3538,0.698497,0.550536,0.134,0.4967,0.1506,0.0677,0.3502,0.7211,0.3161,0.0161,0.6155,0.5647,0.8583"); // Friedrich-ternary-oxide (ICSD #16298)
      }
      // ---------------------------------------------------------------------------
      //DX20201019 [OBSOLETE - moved up ] if(anrl_label=="A2B_oP12_62_2c_c"){
      //DX20201019 [OBSOLETE - moved up ]   vparameters.push_back("3.875,1.64232258065,1.89496774194,0.004,0.758,0.24,0.07,0.24,0.39");
      //DX20201019 [OBSOLETE - moved up ] }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP16_62_cd_c"){
        vparameters.push_back("6.5982,1.11416750023,0.727789397108,0.011,0.415,0.369,0.555,0.174,0.053,0.856");
        vparameters.push_back("4.955,1.4377397,0.90938446,0.05881,0.18533,0.38201,0.93464,0.82626,0.57653,0.31483");  // 002, binary metal-nitride prototype (ICSD #260758)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oP24_62_c_d_cd"){
        vparameters.push_back("6.231,1.78414379714,1.03787514043,0.123,0.0823,0.2579,0.4127,0.1366,0.087,0.5853,0.267,0.0846,-0.088");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP16_62_c_3c"){
        vparameters.push_back("13.855,0.266791771923,0.28601948755,0.398,0.425,0.074,-0.026,0.414,-0.066,0.276,0.49");
        vparameters.push_back("13.52,0.27366864,0.29222633,0.102,0.915,0.437,0.504,0.087,0.479,0.223,0.965");  // 002, binary metal-oxide prototype (ICSD #644065)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP24_62_c_2cd_c"){
        vparameters.push_back("8.884,0.614362899595,0.805155335434,0.8154,0.3419,0.5878,0.6062,0.3192,0.5515,0.437,0.6914,0.4186,0.4702,0.819");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_oC24_63_e_c_cg"){
        vparameters.push_back("9.049,1.21770361366,0.600176815118,0.6699,0.1191,0.8502,0.2174,0.3859");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h"){
        vparameters.push_back("10.11895,1.73974572461,4.16742843872,0.38998,0.0432,0.5987,0.1275,0.13992,0.17412,0.7278,0.20648,0.19723,0.62632,0.1388,0.02884,0.39065,0.0988,0.06788,0.55569,0.59876,0.1328,0.28079,0.54761,0.0678,0.6912,0.07797,0.09224,0.45337,0.1675,0.43146,0.0091,0.18607,0.20422,0.3338,0.3774,0.18828,0.32898,0.174,0.19344,0.2014,0.10123,0.30731,0.03499,0.20655,0.30495,0.49596,0.12681,0.19578,0.33018,0.026,0.31697,0.03483,0.04522,0.2812,0.17213,0.16828,0.2807,0.35936,0.09067");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_oC28_63_efg_c"){
        vparameters.push_back("7.5551,0.860266574896,1.17435904224,0.45686,0.32602,0.13917,0.10039,0.31768,0.28622");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oC20_63_a_cf_c"){
        vparameters.push_back("2.456,3.27442996743,2.48086319218,0.077,0.747,0.631,-0.064"); //DX20210428 - equivalent to V3AsC (http://aflow.org/prototype-encyclopedia/ABC3_oC20_63_c_b_cf.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_63_a_fg_c"){
        vparameters.push_back("5.182,1.52315708221,1.25549980702,0.37,0.25,0.06,0.25,0.47");
        vparameters.push_back("5.765,1.4816999,1.1434519,0.63742,0.7497,0.0434,0.2407,0.0213");  // 002, ternary metal-oxide prototype (ICSD #10431)
        vparameters.push_back("5.568,1.4741379,1.0734555,0.645,0.241,0.974,0.751,0.025");  // 003, ternary metal-oxide prototype (ICSD #27508)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_63_c_fg_c"){
        vparameters.push_back("6.995,0.892780557541,0.999714081487,0.6524,0.15556,0.7025,-0.0819,0.1699,0.0162");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC24_64_2f_f"){
        vparameters.push_back("2.9986,1.44314013206,2.534682852,0.372,0.27,0.6,0.092,0.371,0.615");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oC28_66_l_kl_a"){
        vparameters.push_back("6.2700590867,1.72567783094,1.73046251993,0.833,0.005,0.268,0.737,0.42");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC64_66_gi2lm_2l"){
        vparameters.push_back("8.157,1.00294225818,0.59212945936,0.54552,0.3287,0.39296,0.16145,0.33837,0.89039,0.24136,0.07877,0.42337,0.74171,0.33285,0.66798,0.24948");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC64_66_kl2m_bdl"){
        vparameters.push_back("8.735,2.32364052662,1.67842014883,0.1826,-0.0318,0.1994,0.327,0.1716,0.2894,0.451,0.1302,0.1133,0.3773,0.3708");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_67_ag_b_g"){
        vparameters.push_back("5.0644067238,1.60320657112,1.02379259962,0.3305,0.8198");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_67_b_g_ag"){
        vparameters.push_back("5.2729860526,1.00606865161,1.82912952778,0.765,0.3403");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_67_a_g"){
        vparameters.push_back("5.32495,0.997010300566,1.02896740814,0.26686");
        vparameters.push_back("5.6124272786,0.999376380873,0.8895303257,0.7642");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_oC20_68_a_i"){
        vparameters.push_back("6.44222,1.77662048176,0.991771470083,0.3313,0.1233,0.0801");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oF48_70_f_fg"){
        vparameters.push_back("4.2082,3.45504015969,1.7326647973,0.2495,0.54337,0.29445"); //DX20210428 - this D1_{f} structure is equivalent to Mg2Cu (C_{b}, http://aflow.org/prototype-encyclopedia/AB2_oF48_70_g_fg.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oI14_71_gh_cg"){
        vparameters.push_back("3.29,4.25531914894,0.951367781155,0.375,0.18,0.444");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oI12_71_h_j_g"){
        vparameters.push_back("3.438,3.4554973822,1.37434554974,0.212,0.1232,0.235");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD3_oI48_73_d_e_e_ef"){
        vparameters.push_back("5.7750411261,1.02926406926,3.53056277057,0.63427,0.3734,0.18032,0.311,0.1349,0.6124,0.0967");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oI12_74_h_e"){
        vparameters.push_back("5.158858099,1.56976744188,1.70096899227,-0.047,0.56,0.663"); //DX20210428 - equivalent to CeCu2 (http://aflow.org/prototype-encyclopedia/AB2_oI12_74_e_h.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oI20_74_beh_e"){
        vparameters.push_back("4.39,1.42369020501,3.12528473804,-0.111,0.111,-0.033,0.314");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C12D4_tP76_75_2a2b_2d_12d_4d"){
        vparameters.push_back("9.8880614494,0.922431229769,0.3333,0.0,0.5003,0.1666,0.167,0.348,0.1666,0.333,0.152,0.0,0.208,0.152,0.8333,0.458,0.152,0.8333,0.292,0.348,0.6666,0.042,0.348,0.6666,0.208,0.152,0.5,0.458,0.152,0.5,0.292,0.348,0.3333,0.042,0.348,0.3333,0.208,0.152,0.1666,0.458,0.152,0.1666,0.292,0.348,0.0,0.042,0.348,0.0,0.167,0.348,0.8333,0.333,0.152,0.6666,0.167,0.348,0.5,0.333,0.152,0.3333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP16_76_2a_a_a"){
        vparameters.push_back("3.9810644174,3.85606631503,0.143,0.173,0.1347,0.353,0.141,0.2027,0.3461,0.3475,0.5937,0.1519,0.1578,0.0");
        vparameters.push_back("3.9245,3.8883934,0.145,0.165,0.4605,0.15,0.647,0.637,0.3485,0.3517,0,0.148,0.1533,0.594");  // 002, ternary metal-nitride prototype (ICSD #67327)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_tP40_76_3a_7a"){
        vparameters.push_back("9.046,1.84766747734,0.7435,0.3852,0.0,0.4169,0.733,0.8359,0.026,0.8404,-0.0086,0.79,0.6,0.8105,0.443,0.095,-0.053,0.106,0.473,0.8928,0.357,0.024,0.061,0.629,0.794,0.017,-0.002,0.341,0.7049,0.011,0.29,0.84");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6CD7_tP64_77_2d_6d_d_ab6d"){
        vparameters.push_back("7.6159522234,1.07480314961,0.0,0.035,0.1167,0.1203,0.418,0.392,0.3817,0.132,0.362,0.099,0.309,0.136,0.393,0.161,0.192,0.355,0.442,0.312,0.133,0.03,0.028,0.14,0.14,0.47,0.35,0.32,0.2491,0.7483,0.273,0.2602,0.0188,0.013,0.2316,0.4834,0.184,0.1644,0.2577,0.529,0.3406,0.2365,0.012,0.0182,0.2119,0.275,0.4808,0.3008,0.268");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP48_77_8d_4d"){
        vparameters.push_back("13.47,0.304380103935,0.193,0.14,0.163,0.091,0.059,0.837,0.693,0.14,0.163,0.591,0.059,0.837,0.193,0.64,0.163,0.193,0.64,0.837,0.693,0.64,0.837,0.591,0.559,0.163,0.11,0.14,0.0,0.61,0.14,0.0,0.11,0.64,0.0,0.61,0.64,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_tP88_78_4a_14a_4a"){
        vparameters.push_back("7.1088964505,3.60337042297,0.14561,0.23414,0.69696,0.22057,0.48276,0.79885,0.18673,0.10525,0.29124,0.24906,0.35443,0.19001,0.4586,0.3299,0.02876,0.0936,0.3899,0.1424,0.4187,0.2147,0.16739,0.1085,0.2261,0.23449,0.2999,0.2626,0.32782,0.0777,0.3681,0.7514,0.2725,0.3751,0.65898,0.0556,0.3427,0.02051,0.3251,0.3037,0.52047,0.4811,0.0546,0.59334,0.0026,0.0009,0.31632,0.3813,0.3334,0.82115,0.2817,0.0576,0.71776,0.168,0.0503,-0.08181,0.26434,0.22687,0.42638,0.02701,0.34822,0.57318,0.15408,0.39584,-0.07135,0.37469,0.25769,0.06464");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI20_79_c_2a_c"){
        vparameters.push_back("8.6489709127,0.842525147414,0.0,0.4896,0.337,0.164,0.2196,0.1519,0.1578,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI48_80_2b_4b"){
        vparameters.push_back("9.693,0.617455896007,0.2621,0.5076,0.0299,0.2455,0.4909,0.4804,0.3974,0.1497,0.0077,0.1102,0.3642,-0.0098,0.6086,0.3609,0.5064,0.65,0.1038,0.2484");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP12_81_adg_2h"){
        vparameters.push_back("5.3391020438,1.87980670176,0.25,0.2739,0.234,0.128,0.2289,0.23,0.6373");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI32_82_3g_g"){
        vparameters.push_back("8.954,0.495711413893,0.0775,0.1117,0.2391,0.3649,0.0321,0.9765,0.1689,0.22,0.7524,0.2862,0.0487,0.4807");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_tP10_83_adk_j"){
        vparameters.push_back("6.2840064744,0.638128580522,0.375,0.191,0.109,0.314");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP30_85_ab2g_cg"){
        vparameters.push_back("11.6179902668,0.614391461525,0.6517,0.5428,0.6612,0.5963,0.6531,0.541,0.1258,0.5856,0.1045,0.2524");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP32_86_g_3g"){
        vparameters.push_back("9.997999333,0.499899979999,0.0439,0.20812,0.5354,0.11009,0.22151,0.0295,0.14275,0.66613,0.7153,0.53342,0.06957,0.7593");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tI20_88_f_a"){
        vparameters.push_back("6.407994829,2.01685393259,0.147,0.017,0.298");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI96_88_2f_4f"){
        vparameters.push_back("13.66,0.436603221083,0.1155,0.1249,0.4746,0.1356,0.125,0.0267,-0.0134,0.1262,-0.0046,-0.0251,0.1252,0.5,0.2739,0.1245,-0.0002,0.2631,0.1241,0.5043");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A17BC4D_tP184_89_17p_p_4p_io"){
        vparameters.push_back("16.2503,0.81957871547,0.86772,0.10962,0.1709,0.485,0.7893,0.2334,0.4627,0.7103,0.2146,0.4383,0.6125,0.2885,0.416,0.5631,0.354,0.4242,0.6323,0.3206,0.4579,0.7209,0.234,0.226,0.661,0.3111,0.2275,0.696,0.3148,0.2616,0.794,0.236,0.287,0.818,0.179,0.265,0.733,-0.0696,0.3435,-0.084,-0.078,0.2512,-0.097,-0.004,0.378,0.549,-0.01,0.309,0.593,-0.003,0.238,0.558,-0.006,0.165,0.6,0.2676,0.3463,0.6943,0.1965,0.4869,0.8803,0.096,0.4889,0.7614,0.1143,0.375,0.0177,-0.0146,0.3768,0.8623");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C13D_tP40_90_g_d_cef2g_c"){
        vparameters.push_back("7.3739946979,1.45226471387,-0.0654,0.7769,0.83072,0.6476,0.7846,0.648,0.751,0.001,0.7502,0.5354,0.32757,0.7289,0.8811,0.2635");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C17D4E_tP54_90_a_g_c4g_g_c"){
        vparameters.push_back("9.5601062765,0.748953974895,0.2035,-0.023,0.7639,0.5152,0.407,0.6525,0.8664,0.4465,0.8786,0.838,0.2708,0.6666,0.8964,0.0943,0.6838,0.656,0.2468,0.7207,0.8114,0.2609");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP24_91_d_d_d"){
        vparameters.push_back("3.7620082462,6.71079213192,0.303,0.202,0.019,0.296,0.189,0.08,0.2975,0.1983,0.1795");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB32CD4E8_tP184_93_i_16p_af_2p_4p"){
        vparameters.push_back("11.461,2.92112381119,0.62279,0.369,0.14,0.5403,0.251,0.129,0.6173,0.138,0.647,0.7882,0.215,0.58,0.8645,0.138,0.517,0.591,0.132,0.58,0.5541,0.244,0.605,0.536,0.345,0.554,0.5495,0.343,0.485,0.5842,0.237,0.463,0.6031,0.008,0.366,0.6554,0.077,0.371,0.6892,0.082,0.273,0.7151,0.023,0.17,0.7031,-0.04,0.165,0.6678,-0.05,0.265,0.6408,0.2272,0.1019,0.5636,0.2601,0.5802,0.812,0.1021,0.205,0.5436,0.1936,-0.0676,0.5558,0.404,0.6758,0.8049,0.2809,0.418,0.7912");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A14B3C5_tP44_94_c3g_ad_bg"){
        vparameters.push_back("7.3451315192,1.41592920353,0.8225,-0.0016,0.7086,-0.0024,0.6281,-0.0287,0.6575,0.6242,0.54,0.75,0.5448,0.7312,0.7348,0.7403");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_tP18_94_eg_c_a"){
        vparameters.push_back("4.6863209101,1.96124874637,0.6623,0.7093,0.684,0.707,0.6579");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP24_95_d_d_d"){
        vparameters.push_back("3.7620082462,6.71079213192,0.303,0.202,-0.019,0.296,0.189,-0.08,0.2975,0.1983,0.8205");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8CD_tI24_97_d_k_a_b"){
        vparameters.push_back("5.4068544677,1.92010356944,0.1697,0.3128,0.1237");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_tI44_97_e_2k_cd"){
        vparameters.push_back("9.5317026755,1.33879580768,0.1553,0.0449,0.284,0.3693,0.312,0.1212,0.1191");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI12_98_f_a"){
        vparameters.push_back("7.953376649,0.587954231111,0.44");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C2D_tP26_100_c_abcd_c_a"){
        vparameters.push_back("8.527,0.611047261639,0.7904,0.4646,0.3707,0.32701,0.0,0.1259,0.7949,0.1282,0.4871,0.2924,0.5772,0.3571");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B11C6_tP40_100_ac_bc2d_cd"){
        vparameters.push_back("10.1311200179,0.47843253381,0.0,0.0672,0.18111,0.0154,0.6536,0.6912,0.6174,0.0432,0.5814,0.678,0.6402,0.7288,0.574,0.1742,0.7098,0.5782,0.5334");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B7C2_tP32_101_bde_ade_d"){
        vparameters.push_back("9.8510697809,0.697188102729,0.0,0.4692,0.17136,0.73945,0.2254,0.3402,0.30926,0.0086,0.2384,0.5244,0.2259,0.0352,0.3449,0.0281");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_102_2c_b2c"){
        vparameters.push_back("8.3289849893,0.909833113223,0.5,0.604,0.439,0.623,0.031,0.795,0.725,0.848,0.251");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_tP10_103_a_d"){
        vparameters.push_back("6.5509768136,1.04518394138,0.0,0.144,0.3276,0.242");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B5C4_tP28_104_ac_ac_c"){
        vparameters.push_back("10.6225961282,0.84830508475,0.5,0.8821,0.8116,0.6057,0.3261,0.60942,0.80921,0.00978,0.8116,-0.072,0.1681");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C4_tP22_104_a_2ac_c"){
        vparameters.push_back("9.3940153509,0.981690440703,0.786,0.5,0.0649,0.8297,0.6458,0.286,0.6491,0.8588,0.036");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tP20_105_f_ac_2e"){
        vparameters.push_back("7.7858653925,1.11276650398,0.0,0.0241,0.3351,0.3664,0.1632,0.6068,0.3453,0.229,0.2558");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3D_tP64_106_3c_c_3c_c"){
        vparameters.push_back("10.8312004139,0.492105992062,0.836,0.614,0.2,0.531,0.845,0.11,0.858,0.825,0.08,0.5929,0.62036,0.0343,0.8146,0.6194,0.0688,0.5681,0.8298,0.1322,0.8669,-0.0994,0.0636,0.69285,-0.05677,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B7_tI24_107_ac_abd"){
        vparameters.push_back("7.6400197048,0.760471204184,0.0,0.056,0.04,0.22,0.0,0.243,0.29");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_107_a_a"){
        vparameters.push_back("3.5440505103,1.57477426639,0.0,0.427");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_tI32_108_ac_a2c"){
        vparameters.push_back("8.0549870847,1.94761018001,0.007,0.75,0.109,0.257,0.676,0.114,0.676,0.4");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI12_109_a_a_a"){
        vparameters.push_back("4.2490694941,3.42174629325,0.081,0.666,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI8_109_a_a"){
        vparameters.push_back("3.4517145504,3.38383984705,0.416,0.0"); //DX 20181220 - changed from "0.5416,0.5" to "0.416,0.0"
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC8_tI176_110_2b_b_8b"){
        vparameters.push_back("13.62003,0.668135092213,0.5491,0.8362,0.259,0.8065,0.6386,0.3704,0.1998,0.0869,0.0,0.5687,-0.0955,0.3098,0.6012,0.7807,0.2082,0.5033,0.7931,0.3482,0.6413,0.5083,0.4159,0.8436,0.6054,0.2752,0.8374,0.7122,0.3904,0.7286,0.645,0.3515,0.815,0.5899,0.4761");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP12_111_2n_adf"){
        vparameters.push_back("5.1219931862,1.02616165562,0.205,0.28,0.301,0.622");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP8_111_n_n"){
        vparameters.push_back("4.1305540686,0.998959140196,0.2522,0.7473,0.24415,0.24404");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP12_112_b_n_e"){
        vparameters.push_back("5.4410982776,1.85862633019,0.2334,0.2761,0.6295");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC7D2_tP24_113_e_a_cef_e"){
        vparameters.push_back("7.8338,0.639306594501,0.8201,0.8324,0.4935,0.6407,0.7471,0.6396,0.0642,0.0798,0.1862,0.7856");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP32_114_3e_e"){
        vparameters.push_back("9.6362543341,0.547945205485,0.6,0.743,0.809,0.836,0.555,0.604,0.881,0.899,0.844,0.5125,0.7273,0.563");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tP10_114_e_a"){
        vparameters.push_back("5.2323591487,1.07923706139,0.626,0.768,0.846");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP5_115_g_ag"){
        vparameters.push_back("3.3269188443,1.84881274424,0.253,0.6308");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP12_115_j_egi"){
        vparameters.push_back("8.7863400452,0.701854022742,0.0288,0.0315,0.26398,0.24918,0.24945");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_116_bci_fj"){
        vparameters.push_back("6.1720115185,1.606448477,0.177,0.625,0.655,0.216,0.582");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_117_i_adgh"){
        vparameters.push_back("7.7289660931,0.727131582349,0.73,0.73,0.75,0.52,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP16_118_ei_f"){
        vparameters.push_back("6.9983025398,1.03510852635,0.237,0.15,0.343,0.149,0.509");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_tP32_118_g2i_aceh"){
        vparameters.push_back("5.8229835854,2.43860552978,0.6709,0.675,0.5861,0.73,0.85,0.5515,0.84,0.8,0.15");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI24_119_b2i_af"){
        vparameters.push_back("6.315,2.37529691211,0.372,0.2068,0.2229,0.3067,0.3917");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_119_c_a"){
        vparameters.push_back("5.4790101504,0.558496075934");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2_tI28_120_i_d_e"){
        vparameters.push_back("8.8470588481,0.924381146154,0.856,0.6452,0.6575,0.0851");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4D_tP10_123_gh_a_i_d"){
        vparameters.push_back("3.8757,3.38106664603,0.3336,0.1193,0.2246");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP12_124_a_m_c"){
        vparameters.push_back("6.1884808485,0.816448537725,0.162,0.662");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_tP10_124_a_m"){
        vparameters.push_back("6.4989671731,1.05200800123,0.1425,0.3361");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tP10_125_m_a"){
        vparameters.push_back("6.6398746049,0.899096385539,0.425,0.255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tP12_125_a_b_m"){
        vparameters.push_back("6.3759876428,1.30630489336,0.3822,0.2163");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP28_126_cd_e_k"){
        vparameters.push_back("7.4920241479,1.58609183128,0.11808,0.0865,0.5924,0.1254");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_tP20_127_ehj_g"){
        vparameters.push_back("7.256,0.56684123484,0.2,0.31,0.1,0.2,0.04");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_tP18_128_eh_d_a"){
        vparameters.push_back("7.057532571,1.41383170154,0.2523,0.2217,0.2511");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C_tP40_128_egi_h_e"){
        vparameters.push_back("6.336,2.34690656566,0.366,0.2008,0.165,0.278,0.088,0.198,0.42,0.1");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP28_130_f_c_g"){
        vparameters.push_back("8.5103337343,0.683196239716,0.58,0.5815,0.045,0.136,0.597");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_tP32_130_cg_cf"){
        vparameters.push_back("8.465,1.94329592439,0.2271,0.0095,0.1482,0.57997,0.07997,0.10688");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C4D_tP18_132_e_i_o_d"){
        vparameters.push_back("5.6046,2.34700067801,0.2369,0.26316,0.34803");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C_tP16_132_d_io_a"){
        vparameters.push_back("5.4229923801,1.46597824081,0.3,0.2,0.333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP32_133_h_i2j"){
        vparameters.push_back("9.3810096033,0.497068542799,-0.0329,0.8986,0.658,0.0472");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP24_135_gh_h"){
        vparameters.push_back("8.3218,0.607332548247,0.36248,-0.05789,0.17358,0.13396,0.20929");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tP28_135_gh_h_d"){
        vparameters.push_back("8.4909023894,0.697208809336,0.169,0.114,0.386,0.167,0.175");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP40_137_cdf_3g"){
        vparameters.push_back("8.097,1.41410398913,0.0,0.011,0.511,0.533,0.147,0.467,0.864,0.5,0.603");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_137_d_a"){
        vparameters.push_back("3.64008007,1.4478021978,0.565");
        vparameters.push_back("3.5815,1.4431104,0.3006");  // 002, binary metal-oxide prototype (ICSD #93028)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4_tP18_137_g_b_g"){
        vparameters.push_back("5.0593101041,1.39612571654,0.08,0.1,0.503,0.384");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP6_137_a_d"){
        vparameters.push_back("4.3675,2.8551803091,0.389");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tP12_138_bi"){
        vparameters.push_back("3.388,1.77420306966,0.086,0.107");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI8_139_e_e"){
        vparameters.push_back("4.4795,2.43451278044,0.3356,0.119");
        vparameters.push_back("1.0,2.82842712475,0.125,0.625"); // Lederer-47
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_tI32_140_ah_bk"){
        vparameters.push_back("9.64,0.515560165975,0.17,0.074,0.223");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_tI32_140_ah_cl"){
        vparameters.push_back("5.46,1.91575091575,0.625,0.166,0.15");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_tI16_142_f"){
        vparameters.push_back("8.5939,0.420984651904,0.1405");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B14C3_hP21_143_bd_ac4d_d"){
        vparameters.push_back("7.3813879247,0.611841213925,0.0,0.305,0.497,0.453,0.09,0.302,0.087,0.605,0.536,0.253,0.059,0.583,0.141,0.429,0.08,0.444,0.296,0.071,0.2215,0.2757,0.806");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C_hP11_143_bd_2d_a"){
        vparameters.push_back("6.9672687469,0.526119402977,0.0004,-0.0007,0.8181,0.1915,0.4998,0.5475,0.48,0.0003,0.1859,0.799,0.5002");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP12_143_cd_ab2d"){
        vparameters.push_back("6.5003859369,0.944615384626,0.0,0.5,0.25,0.0548,0.2679,0.25,0.33333,0.16667,0.5,0.0,0.5,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_hP15_144_4a_a"){
        vparameters.push_back("6.2151204722,1.25245374096,0.4904,0.2194,0.2268,0.2226,0.4873,0.1142,0.0775,0.0012,0.0,0.6097,0.0014,0.0018,0.3178,0.0008,0.5062");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_144_a_a"){
        vparameters.push_back("4.0452497415,2.30951792334,0.33,0.16,0.197,0.13,0.32,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C3DE7_hP48_145_2a_3a_3a_a_7a"){
        vparameters.push_back("6.7260126433,2.23669342849,0.567,0.317,0.168,0.432,0.752,0.169,0.6161,-0.0009,0.0,-0.0011,0.6284,0.0045,0.3706,0.371,-0.0043,-0.003,0.266,0.01,0.271,0.002,-0.005,0.73,0.734,-0.02,0.001,0.297,0.165,0.664,0.348,0.09,0.668,0.333,0.238,0.353,0.273,0.167,0.321,0.654,-0.092,0.336,0.666,0.761,0.646,-0.079,0.165,0.203,0.202,0.831");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hR5_146_b_a_a"){
        vparameters.push_back("6.8858547087,1.22475238897,0.47,0.0,0.49,-0.002,-0.14399");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR10_146_2a_2a_2b"){
        vparameters.push_back("6.2648898445,3.16041500397,0.3911,0.7256,0.1132,0.0,0.1454,-0.1886,0.474,0.6151,-0.0128,0.3247");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_hR42_148_2f_4f_f"){
        vparameters.push_back("12.4401,0.661634552777,0.20546,-0.56835,0.60961,0.87265,0.10211,-0.72078,0.66221,-0.37042,-0.03997,0.41696,-0.2507,0.0834,0.70548,0.00275,0.03736,0.37614,-0.32772,0.70785,0.53815,-0.23404,-0.05466");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hR18_148_2f_f"){
        vparameters.push_back("13.0471,0.659280606418,0.7407,-0.76315,0.0238,1.14609,0.38403,-0.59217,0.08714,0.06386,0.3263");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_149_acgi_3l"){
        vparameters.push_back("5.1418156677,2.78268310709,0.33333,0.33333,0.0,0.33333,0.421,0.33333,0.33333,0.246,0.0,0.33333,0.088");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_153_3c_2b"){
        vparameters.push_back("6.017,2.87518697025,0.1111,0.4444,0.1111,0.2222,0.09357,0.4444,0.8888,0.09357,0.77778,0.55558,0.09357");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP9_154_bc"){
        vparameters.push_back("6.9082,0.616557134999,0.876,0.23,0.534,0.051");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_156_b2c_3a2bc"){
        vparameters.push_back("4.239807232,4.83608490569,0.0,0.66667,0.33333,0.75,0.16667,0.5,0.08333,0.41667,0.83333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_156_2ab3c_2ab3c"){
        vparameters.push_back("4.2499813346,4.90823529409,0.375,0.70833,0.5,0.83333,0.04167,0.16667,0.45833,0.79167,0.125,0.33333,0.66667,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_156_ab_ab"){
        vparameters.push_back("4.2794836776,1.67515774714,0.636,0.0,0.894,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B6C2_hP13_157_2ac_2c_b"){
        vparameters.push_back("5.939,1.08233709379,0.264,0.736,0.022,0.5,0.522,0.603,0.215,0.397,0.829");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_158_d_a"){
        vparameters.push_back("6.12,0.924509803922,0.0,0.318,0.027,0.237");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP20_159_bc_2c"){
        vparameters.push_back("7.7488577892,0.813266227893,0.0,0.337,0.154,0.021,0.451,0.061,0.287,0.149,0.282,0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hP28_159_ab2c_2c"){
        vparameters.push_back("7.7479899913,0.724961280333,0.0,0.3649,0.0424,0.3891,0.0408,0.3169,0.3198,0.2712,0.0821,0.5089,0.3172,0.1712,0.2563,0.0274");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C7D_hP26_159_b_ac_a2c_b"){
        vparameters.push_back("6.2653109789,1.6324735851,0.0,0.3057,0.0613,0.4379,0.3425,0.1575,0.2472,0.0015,0.5149,0.3102,0.3347,0.1157,0.0618");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR4_160_b_a"){
        vparameters.push_back("4.405,0.610442678774,0.0,0.52073,-0.02217");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B5_hR26_160_a3bc_a3b"){
        vparameters.push_back("12.71838,0.624030733474,0.672,0.194,0.654,0.01201,0.349,0.58199,0.722,0.35601,1.003,-0.20599,0.998,-0.66,0.355,1.00601,1.033,-0.339,0.28799");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hR3_160_a_a_a"){
        vparameters.push_back("6.1703,0.949908432329,0.0,0.79356,0.25763");
        vparameters.push_back("5.9959433,0.87708454,0,0.408,0.647");  // 002, ternary metal-carbo-nitride prototype (ICSD #85783)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR10_160_5a_5a"){
        vparameters.push_back("3.09,12.2653721683,0.0,0.13333,0.4,0.6,0.86667,0.05,0.18333,0.45,0.65,-0.08333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP5_164_d_ad"){
        vparameters.push_back("3.9381,1.55813717275,0.2467,0.647");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP9_164_bd_c2d"){
        vparameters.push_back("2.89,7.90657439446,0.0607,0.154,0.27263,0.39403");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP4_164_a_b_d"){
        vparameters.push_back("4.0482,1.26777333136,0.271");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_165_bdg_f"){
        vparameters.push_back("7.07,1.00919377652,0.17,0.38,0.69,0.07,0.08");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hR7_166_2c_ac"){
        vparameters.push_back("3.335,7.48635682159,0.29422,0.12967,0.2168");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hR6_166_c_c_c"){
        vparameters.push_back("3.8548,7.95397945419,0.1159,0.3017,0.3815");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR10_167_b_e_a"){
        vparameters.push_back("5.4577,2.40134122433,0.6946");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR24_167_e_e_2e"){
        vparameters.push_back("12.76,0.575235109718,1.1389,0.8113,1.0343,0.3584");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B13C4_hP57_168_d_c6d_2d"){
        vparameters.push_back("15.9361426085,0.244226907628,0.0,0.4965,0.179,0.555,0.035,0.187,0.007,0.349,0.034,0.018,0.168,0.388,0.01,0.195,0.569,0.02,0.093,0.454,0.539,0.273,0.095,0.555,0.2634,0.09,0.086,0.0789,0.4368,0.071");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP72_168_2d_8d_2d"){
        vparameters.push_back("13.759939558,0.609738372094,0.4498,0.113,0.6253,0.1278,0.466,0.1246,0.4154,0.2053,0.0456,0.2066,0.4189,0.5576,0.4218,0.0894,0.8248,0.1514,0.4907,0.3249,0.3746,0.0106,0.0948,0.0083,0.3667,0.4879,0.5693,0.1616,0.0357,0.1505,0.5625,0.5943,0.4451,0.1173,0.0,0.1291,0.4589,0.4993");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP30_169_2a_3a"){
        vparameters.push_back("6.4300116544,2.78071539658,0.013,0.3579,0.12736,0.3339,0.3226,0.29886,0.3347,0.0,0.004,0.0119,0.3343,0.0,0.338,0.0064,0.33823");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP30_170_2a_3a"){
        vparameters.push_back("6.4300116544,2.78071539658,0.013,0.3579,0.87264,0.3339,0.3226,0.70114,0.3347,0.0,-0.004,0.0119,0.3343,0.0,0.338,0.0064,0.66177");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B2C_hP39_171_5c_c_a"){
        vparameters.push_back("6.3199906634,3.05221518986,0.0,-0.004,0.253,0.23633,0.015,0.248,0.11034,0.681,0.192,0.16733,0.448,0.188,0.29533,0.449,0.25,0.048,0.373,0.065,0.50467");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B2C_hP39_172_5c_c_a"){
        vparameters.push_back("6.3199906634,3.05221518986,0.0,-0.004,0.253,0.76367,0.015,0.248,0.88966,0.681,0.192,0.83267,0.448,0.188,0.70467,0.449,0.25,-0.048,0.373,0.065,0.49533");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_173_c_b"){
        vparameters.push_back("7.1329719992,1.03939436422,0.0,0.0337,0.3475,0.146");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hP14_173_bc_c"){
        vparameters.push_back("7.603038022,0.382612126795,0.0,0.3284,0.0313,0.05,0.2314,0.4063,0.013");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B7C2_hP21_174_2j2k_ajk_cf"){
        vparameters.push_back("9.0004021308,0.399102242177,0.4309,0.3719,0.1189,0.2772,0.4163,0.1204,0.0495,0.4359,0.2232,0.124,0.2889,0.4096");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_174_cj_fk_aj"){
        vparameters.push_back("10.7303215747,0.395153774459,0.30167,0.15433,0.03467,0.51733,0.14967,0.31433");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B7C6_hP21_175_ck_aj_k"){
        vparameters.push_back("9.5057625379,0.329093950194,0.36335,0.08627,0.0661,0.221,0.15405,0.51373");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP36_175_jk_jk_jk"){
        vparameters.push_back("11.5799622371,0.317895264089,0.255,0.0598,0.1323,0.416,0.3483,0.0762,0.2334,0.5727,0.4279,0.0715,0.1367,0.5046");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_176_h_bc"){
        vparameters.push_back("7.8700439404,0.50025412961,0.3847,0.0915");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_hP14_176_h_h_c"){
        vparameters.push_back("9.3500327107,0.451336898402,0.1493,0.1701,0.357,0.0462");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_176_h_c"){
        vparameters.push_back("7.4429335392,0.580545478976,0.375,0.083");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP36_177_j2lm_n"){
        vparameters.push_back("12.7835,0.291064262526,0.61855,0.39242,0.79257,0.44445,0.52169,0.86952,0.16458");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_178_b_ac"){
        vparameters.push_back("5.14898393,3.15789473684,0.8361,0.7601,0.5338,0.3099,-0.0053");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_hP6_178_a"){
        vparameters.push_back("2.355,4.43566878981,0.461");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_179_b_ac"){
        vparameters.push_back("5.14898393,3.15789473684,0.8361,0.7601,0.5338,0.3099,0.0053");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP9_181_j_c"){
        vparameters.push_back("4.9977,1.09252256038,0.2072");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP3_183_a_a_a"){
        vparameters.push_back("3.3908401495,1.49568037743,0.608,0.0,0.226");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_183_c_ab"){
        vparameters.push_back("5.3175214551,0.83247442072,0.0,0.513,0.01");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP72_184_d_4d_d"){
        vparameters.push_back("13.7178848276,0.616168537688,0.45652,0.12053,0.0,0.42,0.2069,0.071,0.1224,0.4519,0.2982,0.3629,0.0019,0.0649,0.1538,0.5737,0.0602,0.12298,0.4525,0.12746");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP30_185_cd_c_ab"){
        vparameters.push_back("11.795090915,0.502416278086,0.0,0.377,0.1598,0.2396,0.6647,0.1706,0.1732,0.5056,0.1148");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_185_ab2c_c"){
        vparameters.push_back("6.9593,1.02639633296,0.3213,0.1998,0.2806,0.0765,0.3761,0.4246,0.3322,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_185_c_a"){
        vparameters.push_back("6.1197165237,0.924509803925,0.0,0.305,0.245");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP24_185_c_ab2c"){
        vparameters.push_back("8.7838,1.02449964708,0.2684,0.2311,0.3321,0.25,0.3153,0.5863,0.3518,-0.0769");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_hP20_186_c_b2c"){
        vparameters.push_back("9.85,0.624365482234,0.06,0.815,0.31,0.126,0.25,0.544,0.31");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP4_187_e_fh"){
        vparameters.push_back("2.8065,2.53411722786,0.198");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP10_188_k_c_a"){
        vparameters.push_back("7.2864263258,0.928891999832,0.0041,0.32685");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB9C4_hP28_188_e_kl_ak"){
        vparameters.push_back("6.4953629976,1.43896355826,0.07103,0.48306,0.12023,0.75436,0.22923,0.00127,0.6032");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8BC3D6_hP18_189_bfh_a_g_i"){
        vparameters.push_back("6.62,1.20241691843,0.403,0.444,0.231,0.75,0.222");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9BC3D5_hP18_189_fi_a_g_bh"){
        vparameters.push_back("6.6,1.19696969697,0.378,0.43,0.266,0.755,0.236");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP18_190_gh_bf"){
        vparameters.push_back("7.946,0.789076264787,0.0225,0.294,0.612,0.01");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_hP16_190_bdh_g"){
        vparameters.push_back("6.9236970109,1.22634969237,0.3313,0.0628,0.6682");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP24_190_i_afh"){
        vparameters.push_back("5.9699820408,1.96984924623,0.52,0.6683,0.6653,0.3786,0.3233,0.623");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C18D6_hP58_192_c_f_lm_l"){
        vparameters.push_back("9.214,0.997829390059,0.3103,0.2369,0.3876,0.1159,0.4985,0.1456,0.1453");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP72_192_m_j2kl"){
        vparameters.push_back("13.7705812079,0.608458538779,0.6373,0.7895,0.4221,0.6688,0.5445,0.1221,0.4551,0.686");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_hP16_193_dg_g"){
        vparameters.push_back("6.9104160691,0.696671490596,0.2358,0.5992");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP16_194_gh_ac"){
        vparameters.push_back("5.096,1.6295133438,-0.16667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hP28_194_ahk_ch"){
        vparameters.push_back("7.656,0.991771159875,0.533,0.872,0.196,0.439");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B3C_hP26_194_hk_h_a"){
        vparameters.push_back("7.513,1.03087980833,0.458,0.12,0.201,-0.067");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12BC4_cP34_195_2j_ab_2e"){
        vparameters.push_back("8.0357772599,0.2493,0.7507,0.14225,0.5,0.35579,0.0002,0.14286,0.35831");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B2C_cF60_196_h_bc_a"){
        vparameters.push_back("9.9799933334,0.0,0.0625,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_cP20_198_a_a_b"){
        vparameters.push_back("6.57,0.417,0.064,0.303,0.592,0.5");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B11_cP39_200_f_aghij"){
        vparameters.push_back("8.5520223662,0.18,0.34,0.265,0.278,0.157,0.257");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_cP60_201_be_fh_g"){
        vparameters.push_back("9.5599167841,0.3889,0.6111,0.0972,0.0389,0.0972,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6C_cF104_202_h_h_c"){
        vparameters.push_back("10.6100296668,0.5827,0.6359,0.638,0.72");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cF240_202_h2i"){
        vparameters.push_back("14.26,0.249,0.052,0.105,0.085,0.22,0.185,0.052,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD3E6_cF208_203_e_c_d_f_g"){
        vparameters.push_back("13.9898,0.2826,-0.099,0.2257,0.2665,0.3531");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C6D16E_cF232_203_e_d_f_eg_a"){
        vparameters.push_back("13.9038,0.28207,0.06362,0.34379,0.26626,0.22529,0.35333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C16_cF160_203_a_bc_eg"){
        vparameters.push_back("16.6600284675,0.20516,0.01201,0.111,0.42978");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C6_cP264_205_2d_ab2c2d_6d"){
        vparameters.push_back("15.263,0.2561,0.375,0.2526,0.0133,0.0197,0.2444,0.2335,0.0046,0.1386,0.3763,0.1272,0.38,0.3838,0.1209,0.2777,0.1241,0.0103,0.4835,0.1315,0.2536,0.2664,0.2841,0.1049,0.235,0.4047,0.2921,0.3491,-0.0385,-0.1074,0.1509,-0.0104,-0.0242");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_cP240_205_10d"){
        vparameters.push_back("14.04078,0.2294,-0.0325,0.101,0.2467,-0.054,0.0061,0.2081,0.0646,0.1289,0.2066,0.8599,-0.036,0.171,-0.0963,0.159,0.2236,0.1122,-0.0371,0.2439,0.0192,-0.0636,0.2053,0.1349,0.0616,0.1503,0.7983,0.0202,0.1323,0.8207,0.1186");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_cI96_206_c_e_ad"){
        vparameters.push_back("9.46,0.115,0.205,0.16,0.382,0.11");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A17B15_cP64_207_acfk_eij"){
        vparameters.push_back("10.6058825779,0.2422,0.2622,0.3319,0.2701,0.142,0.1539,0.3498");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_cP16_208_j_b"){
        vparameters.push_back("6.31,0.184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2CD6E_cP64_208_m_ad_b_m_c"){
        vparameters.push_back("10.3701312618,0.057,0.25,0.23,0.235,0.25,0.458");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A24BC_cF104_209_j_a_b"){
        vparameters.push_back("7.7099775082,0.043,0.109,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B36CD12_cF488_210_h_3h_a_fg"){
        vparameters.push_back("16.4321054599,0.12536,0.51382,0.55423,0.58845,0.50025,0.5008,0.59,0.3557,0.5123,0.5411,0.1708,0.5347,0.2125,0.5888");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B6C_cF608_210_4h_2h_e"){
        vparameters.push_back("16.4321054599,0.3771,0.0157,0.2009,0.1224,0.5287,0.1425,0.0068,0.0949,0.2596,0.2225,0.6117,0.1785,0.0224,0.0928,0.2607,0.1616,0.0194,0.0793,0.3503");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cI72_211_hi_i"){
        vparameters.push_back("9.68882,0.37338,0.15866,0.38235");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cP12_212_c_a"){
        vparameters.push_back("6.54,0.428");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_cI56_214_g_h_a"){
        vparameters.push_back("12.31504,0.108,0.384");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_cI48_214_f_a_e"){
        vparameters.push_back("10.38,0.266,0.365");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B9_cP52_215_ei_3efgi"){
        vparameters.push_back("8.7068,0.1157,0.8296,0.3253,0.6066,0.3534,0.8549,0.8113,0.5332,0.3153,0.0322");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_cF16_216_c_d_b_a"){
        vparameters.push_back("6.465");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_cP16_218_c_e_a"){
        vparameters.push_back("6.0258147002,0.6486");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7BC3D13_cF192_219_de_b_c_ah"){
        vparameters.push_back("12.0986,0.0808,0.0987,0.0214,0.1821");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A15B4_cI76_220_ae_c"){
        vparameters.push_back("9.718,-0.042,0.12,0.16,-0.04");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cI28_220_c_a"){
        vparameters.push_back("8.6,0.08333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C6_cP33_221_cd_ag_fh"){
        vparameters.push_back("7.624,0.3,0.2,0.2");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C16_cP96_222_ce_d_fi"){
        vparameters.push_back("11.1199004528,0.0,0.125,0.084,0.166,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A23B6_cF116_225_bd2f_e"){
        vparameters.push_back("12.523,0.203,0.178,0.378");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_cF36_225_e_c_a"){
        vparameters.push_back("9.725,0.24");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB13_cF112_226_a_bi"){
        vparameters.push_back("12.2836,0.1806,0.1192");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C7_cF88_227_c_d_af"){
        vparameters.push_back("10.2663,0.4157");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_cF56_227_ad_e"){
        vparameters.push_back("8.0835,0.2642");
        vparameters.push_back("8.085,0.7637");  // 002, binary metal-oxide prototype (ICSD #63164)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BCD6_cF416_228_eg_c_b_h"){
        vparameters.push_back("22.2649775014,0.19,0.425,0.05,0.18,0.28");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_cF224_228_h_c"){
        vparameters.push_back("15.4800297791,0.043,0.138,0.278");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B10_cI52_229_e_fh"){
        vparameters.push_back("8.9822,0.3538,0.13505,0.3045");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cI10_229_c_a"){
        vparameters.push_back("6.186");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B3_cI40_229_df_e"){
        vparameters.push_back("8.735,0.342,0.156");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C12D3_cI160_230_a_c_h_d"){
        vparameters.push_back("11.4597,0.3471,0.4664,0.0512");
      }
    }
    // ---------------------------------------------------------------------------
    // miscellaneous structures
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="misc"){
      // SQS structures
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_aP16_2_4i_4i"){
        vparameters.push_back("5.4244498076,1.04446593575,1.78376517004,72.976133815,87.0786711125,79.9750121379,0.375,0.75,0.375,0.875,0.0,0.375,0.375,0.5,0.875,0.625,0.75,0.625,0.875,0.5,0.375,0.125,0.75,0.125,0.875,0.75,0.875,0.375,-0.0,0.875");
        vparameters.push_back("6.65,0.83308271,1.0541353,90.1,90.4,90.9,0.065,0.251,0.106,0.668,0.256,0.116,0.061,0.259,0.64,0.605,0.246,0.63,0.94,0.25,0.37,0.5,0.25,0.36,0.25,0.25,0.87,0.8,0.25,0.86");  // 002, binary metal-oxide prototype (ICSD #32561)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B11_mP16_6_2abc_2a3b3c"){
        vparameters.push_back("6.5421326204,1.0,1.0,90.0,0.5,-0.0,0.5,0.5,0.0,-0.0,0.0,0.5,0.0,0.5,0.0,-0.0,0.5,-0.0,0.5,0.5,0.25,0.25,0.25,0.75,0.25,0.25,0.75,0.25,0.75,0.25,0.25,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC32_8_4a_12a"){
        vparameters.push_back("13.8779590179,0.333333333333,0.666666666667,109.471220634,-0.0,0.5,0.5,0.75,0.375,0.9375,0.125,0.8125,0.875,0.1875,0.625,0.0625,0.75,0.375,0.875,0.6875,-0.0,-0.0,0.5,0.25,0.625,0.5625,0.75,0.875,0.25,0.125,0.375,0.4375,0.125,0.3125,0.25,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC32_8_4a_4a4b"){
        vparameters.push_back("9.2519726786,1.0,0.707106781188,90.0,0.0,0.0,0.0,0.5,0.25,0.25,0.5,0.5,0.25,0.75,0.5,0.0,0.75,0.75,0.75,0.25,0.75,0.25,0.0,0.75,0.25,0.5,0.0,0.25,0.75,0.0,0.25,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B13_oC32_38_ac_a2bcdef"){
        vparameters.push_back("6.5421326204,1.41421356237,1.41421356237,0.5,0.0,0.5,0.0,0.75,0.25,0.75,0.75,0.25,0.25,0.25,0.25,0.75,0.75,0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oC32_38_abce_abcdf"){
        vparameters.push_back("6.5421326204,1.41421356237,1.41421356237,0.5,-0.0,-0.0,0.5,0.25,0.75,0.25,0.25,0.75,0.25,0.75,0.25,0.25,0.75,-0.0");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB7_hR16_166_c_c2h"){
        vparameters.push_back("9.2519726786,1.22474487138,0.875,0.625,1.625,0.1250000001,1.125,1.6249999999");
      }
      // kesterite
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD4_tI16_82_ac_b_d_g"){
        vparameters.push_back("5.427,2.00313248572,0.7434,0.256,0.6278");
      }
      // misc structures (Y. Lederer)
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC8_12_a_di"){
        vparameters.push_back("1.0,0.301511344577,0.522232967868,100.024987862,0.75,0.25");
        vparameters.push_back("5.627,0.58983473,0.88484095,72.6,0.6048,0.2603");  // 002, binary metal-nitride prototype (ICSD #34675)
        vparameters.push_back("5.328,0.60041291,0.8871997,77.38,0.3756,0.7456");  // 003, binary metal-nitride prototype (ICSD #181559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_12_i_i"){
        vparameters.push_back("1.0,0.301511344581,0.522232967872,100.024987862,0.625,0.875,0.875,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_mC16_12_a_di_2i"){
        vparameters.push_back("1.0,0.301511344577,0.522232967868,100.024987862,0.75,0.25,0.125,0.375,0.625,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC16_12_i_i_adi"){
        vparameters.push_back("1.0,0.301511344581,0.522232967872,100.024987862,0.625,0.875,0.875,0.625,0.75,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oP4_47_a_ct"){
        vparameters.push_back("1.0,1.41421356239,2.0,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP4_47_cr_a"){
        vparameters.push_back("1.0,1.41421356238,2.82842712476,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP8_47_a_ct_egs"){
        vparameters.push_back("1.0,1.41421356239,2.0,0.25,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_oP8_47_eq_g_bdt"){
        vparameters.push_back("1.0,1.41421356238,2.82842712476,0.75,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_51_e_e"){
        vparameters.push_back("1.0,0.707106781182,2.0,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP8_51_e_e_2f"){
        vparameters.push_back("1.0,0.707106781182,2.0,0.125,0.625,0.375,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_59_a_a"){
        vparameters.push_back("1.0,1.41421356238,2.0,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP8_59_a_a_2b"){
        vparameters.push_back("1.0,1.41421356238,2.0,0.875,0.375,0.875,0.375");
        vparameters.push_back("2.852,1.6725105,2.2124825,0.375,0.875,0.875,0.375");  // 002, ternary metal-oxide prototype (ICSD #16271)
        vparameters.push_back("2.806,1.6215253,2.0481112,0.125,0.643,0.136,0.6");  // 003, ternary metal-oxide prototype (ICSD #84642)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_c_g"){
        vparameters.push_back("1.0,1.41421356236,0.707106781182,0.625,0.125,0.25,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oC12_65_a_i_cj"){
        vparameters.push_back("1.0,2.99999999998,0.707106781176,0.3333333333,0.8333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_oC16_65_ai_b_q"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.75,0.75,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_oC16_65_bj_a_eh"){
        vparameters.push_back("1.0,1.41421356236,0.707106781182,0.75,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oC16_65_a_bf_hi"){
        vparameters.push_back("1.0,1.41421356239,0.5,0.75,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC6_65_a_i"){
        vparameters.push_back("1.0,2.99999999998,0.707106781176,0.3333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC8_65_ai_b"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oC8_65_bj_a"){
        vparameters.push_back("1.0,1.41421356236,0.707106781182,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_65_i_i"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_65_i_i_fh"){
        vparameters.push_back("1.0,1.99999999999,0.499999999994,0.75,0.875,0.375");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oI12_71_a_e_df"){
        vparameters.push_back("1.0,0.471404520797,0.333333333338,0.6666666667,0.6666666667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_a_b"){
        vparameters.push_back("1.0,2.00000000002");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_a_c"){
        vparameters.push_back("1.0,0.707106781182");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP3_123_g_a"){
        vparameters.push_back("1.0,3.00000000002,0.6666666667");
        vparameters.push_back("1.7825,1.4129032,0.638");  // 002, binary metal-oxide prototype (ICSD #92091)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_abc_d"){
        vparameters.push_back("1.0,1.41421356238");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_ag_b"){
        vparameters.push_back("1.0,4.00000000002,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_cf_a"){
        vparameters.push_back("1.0,0.499999999994");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tP4_123_a_bh"){
        vparameters.push_back("1.0,2.82842712477,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_a_b_h"){
        vparameters.push_back("1.0,2.00000000002,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_a_c_e"){
        vparameters.push_back("1.0,0.707106781182");
        vparameters.push_back("3.837,0.95282773");  // 002, metallic ternary prototype (ICSD #42564)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_123_a_d_bc"){
        vparameters.push_back("1.0,1.41421356238");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_123_g_g"){
        vparameters.push_back("1.0,4.00000000002,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_tP6_123_g_b_ch"){
        vparameters.push_back("1.0,3.00000000002,0.1666666667,0.6666666667");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tP8_123_abc_d_i"){
        vparameters.push_back("1.0,1.41421356238,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tP8_123_ag_b_2h"){
        vparameters.push_back("1.0,4.00000000002,0.25,0.625,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tP8_123_cf_a_k"){
        vparameters.push_back("1.0,0.499999999994,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_tP8_123_a_bh_cdg"){
        vparameters.push_back("1.0,2.82842712477,0.25,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_123_h_h_abg"){
        vparameters.push_back("1.0,4.00000000002,0.25,0.625,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_129_c_c_2c"){
        vparameters.push_back("1.0,2.82842712472,0.875,0.375,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tI8_139_ae_b"){
        vparameters.push_back("1.0,2.82842712475,0.25");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_tI12_139_a_e_be"){
        vparameters.push_back("1.0,4.24264068707,0.3333333333,0.8333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_tI16_139_ae_b_g"){
        vparameters.push_back("1.0,2.82842712475,0.25,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_tI16_139_a_bd_ce"){
        vparameters.push_back("1.0,2.0,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_139_e_e_cd"){
        vparameters.push_back("1.0,2.82842712475,0.125,0.625");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_141_a_b_e"){
        vparameters.push_back("1.0,1.99999999997,0.625");
        vparameters.push_back("4.049,2.1590516,0.625");  // 002, ternary metal-oxide prototype (ICSD #31149)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP3_164_d_a"){
        vparameters.push_back("1.0,2.44948974278,0.3333333333");
        vparameters.push_back("1.0,1.22474487139,0.3333333333");
        vparameters.push_back("3.38,1.2810651,0.23");  // 003, binary metal-nitride prototype (ICSD #262746)
        vparameters.push_back("3.08,1.3603896,0.34");  // 004, binary metal-oxide prototype (ICSD #76431)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_164_d_a_bd"){
        vparameters.push_back("1.0,2.44948974278,0.3333333333,0.8333333333");
        vparameters.push_back("4.704,1.6001276,0.742,0.342");  // 002, metallic ternary prototype (ICSD #616769)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_164_d_b_ad"){
        vparameters.push_back("1.0,1.22474487139,0.8333333333,0.3333333333");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR4_166_bc_a"){
        vparameters.push_back("1.0,4.89897948557,0.25");
        vparameters.push_back("3.6459845,4.1752552,0.42402");  // 002, binary metal-nitride prototype (ICSD #1144)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR4_166_a_bc"){
        vparameters.push_back("1.0,9.79795897119,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR4_166_c_c"){
        vparameters.push_back("1.0,9.79795897107,0.625,0.125");
        vparameters.push_back("1.0,4.89897948557,0.375,0.875");
        vparameters.push_back("3.0580228,5.2517983,0.33333333,0.085833333");  // 003, binary metal-boride prototype (ICSD #1)
        vparameters.push_back("2.7849821,6.7813387,0.963,0.771");  // 004, binary metal-carbide prototype (ICSD #188285)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_hR8_166_bc_a_2c"){
        vparameters.push_back("1.0,4.89897948557,0.25,0.625,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_hR8_166_a_bc_2c"){
        vparameters.push_back("1.0,9.79795897119,0.75,0.375,0.125");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR8_166_c_c_abc"){
        vparameters.push_back("1.0,9.79795897107,0.625,0.125,0.75");
        vparameters.push_back("1.0,4.89897948557,0.375,0.875,0.75");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_cP8_221_a_c_bd"){
        vparameters.push_back("1.0");
      }
      // -------------------------------------------------------------------------
      // oxide prototypes (from R. Friedrich)
      // -------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC20_12_b2i_c2i"){
        vparameters.push_back("9.7915972616,0.443468950753,0.626873661672,72.47,0.825,0.343,0.821,0.833,0.66,0.17,0.658,0.669");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_mC20_12_2i_3i"){
        vparameters.push_back("12.8749583858,0.24865727853,0.999198732944,152.527888576,0.7041,0.7946,0.65536,0.31402,0.9453,0.1098,0.3899,0.5632,0.7607,0.2566"); //DX20210427 - equivalent to beta-Ga2O3 (http://aflow.org/prototype-encyclopedia/A2B3_mC20_12_2i_3i.html, part 3)
        vparameters.push_back("11.737,0.24341825,0.96137003,141.8629,0.3177,0.1547,0.698,0.3577,0.5851,0.4465,0.0366,0.8467,0.2674,0.7938");  // 002, binary metal-oxide prototype (ICSD #263119)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_mC32_12_5i_3i"){
        vparameters.push_back("9.7349670426,0.385643337648,1.40884015867,133.94703436,0.0136,0.2506,0.8276,0.8736,0.5586,0.6116,0.6158,0.9368,0.267,0.574,0.442,0.2461,0.9228,0.0529,0.2983,0.4337");
        vparameters.push_back("9.752,0.38986874,0.96821165,88.08,0.676,0.94,0.241,0.755,0.588,0.655,0.953,0.842,0.866,0.559,0.128,0.956,0.7786,0.7331,0.0538,0.6341");  // 002, binary metal-oxide prototype (ICSD #26492)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mP16_14_3e_e"){
        vparameters.push_back("8.3369865333,0.56178500085,0.838076298071,144.831245943,0.958,0.961,0.756,0.442,0.046,0.243,0.255,0.673,0.243,0.2563,0.7319,0.9681");
        vparameters.push_back("6.44,0.93322981,1.3586957,122,0.718,0.721,0.659,0.67,0.749,0.954,0.903,0.891,0.727,0.771,0.243,0.918");  // 002, binary metal-oxide prototype (ICSD #6094)
        vparameters.push_back("6.441,0.93619003,1.3578637,122.25,0.7198,0.7875,0.611,0.6762,0.7521,0.9448,0.8993,0.9329,0.7083,0.7712,0.2434,0.9152");  // 003, binary metal-oxide prototype (ICSD #47164)
        vparameters.push_back("6.72,0.87425595,1.3443452,124.8,0.7137,0.7804,0.5854,0.6844,0.7487,0.921,0.8972,0.931,0.6977,0.7793,0.2472,0.9075");  // 004, binary metal-oxide prototype (ICSD #180568)
        vparameters.push_back("8.2005,0.91823669,0.89043351,152.0345,0,0.97,0.75,0,0.75,0,0.5,0.72,0.25,0.053,0.729,0.809");  // 005, binary metal-oxide prototype (ICSD #647640)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_mP20_14_2e_3e"){
        vparameters.push_back("7.592734201,1.09598618914,0.998035732622,133.736291174,0.524,0.816,0.158,0.04,0.959,0.263,0.778,0.702,0.071,0.235,0.948,0.11,0.27,0.969,0.758");
        vparameters.push_back("4.22,2.26777251185,2.88570616114,109.77475,0.124,0.19,0.267,0.841,0.323,0.484,0.895,0.12,0.395,0.342,0.397,0.355,0.715,0.293,0.125");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_15_a_e"){
        vparameters.push_back("4.6509294963,0.729773531807,1.35357481762,127.166266704,0.333");
        vparameters.push_back("4.6797,0.73325213,1.3606214,127.2378,0.1678");  // 002, binary metal-oxide prototype (ICSD #67850)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_mC20_15_2f_e"){
        vparameters.push_back("7.7087709852,0.481394604967,0.920353982304,63.4,0.7399,0.379,0.462,0.19,0.885,0.524,0.425");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5_oP28_19_2a_5a"){
        vparameters.push_back("4.4859196717,1.82663207955,1.86900129701,0.3784,0.3475,0.4015,0.384,0.9674,0.7836,0.0951,0.3166,0.5403,0.8671,0.6061,0.7453,0.3371,0.4845,0.8324,0.7112,0.2743,0.749,0.55,0.4939,0.5257");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_33_4a_2a"){
        vparameters.push_back("5.6320569384,0.884841795441,2.16335540838,0.6753,0.1608,0.593,0.6454,0.8379,0.906,0.834,0.699,0.6957,0.9256,0.188,0.8093,0.9786,0.9654,0.0,0.6272,0.9994,0.7473"); //DX20210428 - equivalent to alpha-Sb2O4 (cervantite, http://aflow.org/prototype-encyclopedia/A2B_oP24_33_4a_2a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_40_b_a2b"){
        vparameters.push_back("5.743,1.48999,0.833885,0.3841,0.90324,0.5,0.7323,0.3755,0.8922,0.8284");
        vparameters.push_back("5.743,1.54199895525,0.833884729236,0.444,0.097,0.306,0.778,0.0,0.222,0.0");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_oP14_59_a2e_e"){
        vparameters.push_back("3.2077,3.2300714,1.2255822,0.999,0.3957,0.469,0.5689,0.997,0.39882,0.1083"); //DX20210428 - equivalent to O5V2 (shcherbinaite (revised), http://aflow.org/prototype-encyclopedia/A5B2_oP14_59_a2f_f.html, part 3)
        vparameters.push_back("3.571,3.2327079,1.2273873,0,0.5691,0.0028,0.8955,0.5303,0.39877,0.8914");  // 002, binary metal-oxide prototype (ICSD #43132)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC16_64_e_f"){
        vparameters.push_back("8.7597573116,1.03356994972,0.940698018333,0.6554,0.9104,0.0725");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_tP28_135_gh_dh"){
        vparameters.push_back("8.7553254034,0.751693002255,0.33,0.4,0.1,0.14,0.165");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP15_144_2a_3a"){
        vparameters.push_back("3.9844538913,1.92300578035,0.15,0.61,0.02,0.77,0.18,0.26,0.15,0.95,0.0,0.79,0.33,0.07,0.23,0.72,0.56");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hR3_166_c_a"){
        vparameters.push_back("5.73015,4.46193,0.744");
        vparameters.push_back("3.6378513,5.1624334,0.268");  // 002, binary metal-nitride prototype (ICSD #22231)
        vparameters.push_back("4.0290164,5.5658926,0.26647");  // 003, binary metal-nitride prototype (ICSD #409851)
        vparameters.push_back("3.6240359,2.1992272,0.25");  // 004, binary metal-oxide prototype (ICSD #108886)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hR6_166_c_2c"){
        vparameters.push_back("3.5572963728,10.7622298067,0.083,0.377,0.21");
        vparameters.push_back("3.8548,7.95397945419,0.083,0.183,0.35");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_189_fg_eh"){
        vparameters.push_back("7.5487706231,0.71987757732,0.668,0.295,0.632,0.168");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_194_ac_f"){
        vparameters.push_back("3.7639841644,2.4272070374,0.6497");
        vparameters.push_back("2.952,3.8109756,0.625");  // 002, binary metal-nitride prototype (ICSD #76008)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_cI80_199_a2b_2c"){
        vparameters.push_back("10.7844999416,0.5,0.229,0.708,0.625,0.125,0.875,0.625,0.875,0.875");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_cF80_227_f_e"){
        vparameters.push_back("10.7531,0.24,0.895");
        vparameters.push_back("11.1519,0.26027,0.43875");  // 002, (part 3)
      }
      // -------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B4C_aP18_2_4i_4i_i"){
        vparameters.push_back("6.63976352,1.00001361035,1.50518591105,71.844867398,80.9249357443,67.4411683666,0.7435,0.7649,0.0197,0.2323,0.4828,0.5617,0.2681,0.0455,0.3887,0.236,0.5962,0.1682,0.9572,0.3436,0.1926,0.9557,0.8502,0.3514,0.6442,0.1527,0.0884,0.5375,0.2907,0.3741,0.7148,0.2162,0.2479");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_aP22_2_2i_7i_2i"){
        vparameters.push_back("7.3384982604,0.986316404633,0.950710549193,96.331,116.158,86.387,0.82535,0.34153,0.25581,0.23479,0.93732,0.27034,0.5623,0.6988,0.9539,0.6304,0.3997,0.4639,0.6025,0.9461,0.341,0.9238,0.6849,0.3595,0.8578,0.0087,0.1317,0.173,0.2822,0.2504,0.232,0.6233,0.107,0.34536,0.44585,0.27538,0.73214,0.83509,0.19625");
        vparameters.push_back("10.4825,0.46040544,0.50745528,81.42,79.6425,102.8254,0.9163,0.4905,0.7907,0.58313,0.5771,0.21504,0.2925,0.1055,0.6186,0.0477,0.7911,0.917,0.2349,0.8135,0.2409,0.0666,0.3044,0.5772,0.5412,0.7648,0.8863,0.7877,0.6364,0.0493,0.5614,0.2891,0.553,0.36066,0.96292,0.34888,0.1386,0.1407,0.786");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_aP30_2_3i_9i_3i"){
        vparameters.push_back("8.440104917,0.923591309392,0.891430518055,90.055,95.217,103.426,0.19831,0.42266,0.7606,0.20241,0.92919,0.76401,0.50333,0.7504,0.52691,0.3034,0.4616,0.4628,0.3014,0.9385,0.4641,0.5705,0.7688,0.1988,0.9832,0.3739,0.2655,0.9819,0.8677,0.2648,0.4018,0.7266,0.8296,0.2183,0.1785,0.2254,0.2713,0.8704,0.0938,0.2735,0.5126,0.0931,0.1851,0.3875,0.2684,0.1849,0.9542,0.2691,0.3973,0.7236,0.0561"); //DX20210428 - equivalent to CaSiO3 (wollastonite, http://aflow.org/prototype-encyclopedia/AB3C_aP30_2_3i_9i_3i.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_aP32_2_4i_10i_2i"){
        vparameters.push_back("7.3087649311,1.14261433201,0.821778898093,90.89,101.58,105.99,0.3463,0.7028,0.4634,0.2875,0.6959,0.9481,0.0961,0.3757,0.6278,0.1278,0.9318,0.1815,0.1019,0.1479,0.1134,0.1351,0.6991,0.197,0.2641,0.4568,0.9447,0.2832,0.9314,0.9431,0.1337,0.6232,0.6318,0.2765,0.4421,0.4144,0.2948,0.9435,0.4818,0.4953,0.2777,0.2314,0.1057,0.1465,0.6744,0.5003,0.2278,0.7656,0.2959,0.0562,0.72,0.2808,0.3359,0.1702"); //DX20210427 - equivalent to Al2O5Si (kyanite, http://aflow.org/prototype-encyclopedia/A2B5C_aP32_2_4i_10i_2i.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mP28_4_4a_8a_2a"){
        vparameters.push_back("9.0939150001,1.04504737499,0.608139562228,92.875,0.807,0.829,0.272,0.204,0.832,0.263,0.302,0.661,0.774,0.321,0.171,0.214,0.7387,0.1707,0.568,0.2717,0.3202,0.4218,0.6785,0.4914,0.642,0.7422,0.9861,0.0989,0.8129,0.3034,0.0693,0.7893,0.6712,0.078,0.5051,0.2076,0.129,0.0106,0.8783,0.3432,0.5032,-0.0,0.747,0.975,0.9939,0.7889");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_mC60_5_ab8c_ab2c_3c"){
        vparameters.push_back("13.2962734468,0.579323216751,0.88233384728,68.42,0.126,0.583,0.395,0.925,0.458,0.177,0.139,0.393,0.968,0.991,0.37,0.658,0.858,0.777,0.165,0.139,0.451,0.315,0.635,0.89,0.089,0.511,0.381,0.841,0.35,0.287,0.857,0.63,0.384,0.71,0.999,0.37,0.798,0.498,0.9127,0.0,0.2494,0.2545,0.022,0.2458,0.5872,0.0207,0.2516");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5C_mC54_8_3a3b_9a3b_3a"){
        vparameters.push_back("13.993228409,0.578095627303,0.95071704775,134.2304365,0.0445,0.3334,0.9668,0.6832,0.5007,0.0005,0.2994,0.4493,0.9984,0.4969,0.6404,0.2828,0.8665,0.803,0.7033,0.5434,0.7612,0.1446,0.2517,0.8741,0.3501,0.7145,0.9192,0.057,0.0021,0.9993,0.5651,0.3477,0.4285,0.6427,0.1992,0.2336,0.6793,0.3112,0.2665,0.3313,0.7507,0.2478,0.0016,0.4657,0.8234,0.2857,0.5968,0.6856,0.0533,0.035,0.3071,0.7119");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C3_mP24_11_2e_7e_3e"){
        vparameters.push_back("9.7906148128,0.416420361249,0.938259441704,101.57,0.595,0.682,0.154,0.508,0.195,0.221,0.473,0.14,0.645,0.438,0.885,0.314,0.745,0.997,0.313,0.791,0.031,0.905,0.2806,0.0278,0.673,0.2467,0.9811,0.142");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2_mC18_12_a_3i_i"){
        vparameters.push_back("9.7742615012,0.36510934394,0.69960238568,75.2,0.6434,0.8818,0.7535,0.3934,0.9126,0.7219,0.7658,0.6691");
        vparameters.push_back("9.1042,0.37740823,0.9940467,137.3359,0.4616,0.1133,0.3157,0.2764,0.8722,0.5652,0.53846,0.3436");
        vparameters.push_back("9.2545,0.37866984,0.99218758,137.9112,0.1881,0.2237,0.0431,0.39,0.3729,0.0653,0.961,0.155");  // 003, ternary metal-oxide prototype (ICSD #188911)
        vparameters.push_back("12.46083,0.270848731585,0.836461134611,116.31902,-0.00589,0.76131,0.65848,0.53954,0.21851,0.12521,0.67999,0.28976");  // 004, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mC48_12_gi_hi_2i3j"){
        vparameters.push_back("9.7995229272,0.904117589797,0.6838314027,73.04,0.3216,0.249,0.2996,0.3569,0.2291,0.9043,0.3587,0.0391,0.2983,0.6449,0.5415,0.8467,0.696,0.3587,0.6439,0.6088,0.6337,0.6552,0.9717");
        vparameters.push_back("8.8789,0.912805,0.986699,133.694,0.7148,0.1906,0.3547,0.6391,0.8447,0.1555,0.09,0.4376,0.5669,0.9156,0.73601,0.8485,0.23951,0.4244,0.3344,0.0772,0.4173,0.8466,0.5815");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mP12_13_f_2g_e"){
        vparameters.push_back("6.912898268,0.837946008451,0.729410860817,136.193216009,0.1825,0.327,0.22,0.894,0.281,0.256,0.622,0.855"); //DX20210428 - equivalent to MgO4W (huanzalaite, H0_{6}, http://aflow.org/prototype-encyclopedia/AB4C_mP12_13_f_2g_e.html, part 3)
        vparameters.push_back("6.5113,0.84124522,0.72976211,136.867,0.82,0.653,0.78,0.11,0.82,0.74,0.38,0.44");
        vparameters.push_back("6.7849,0.83567923,0.72661351,136.2724,0.182,0.32,0.2,0.96,0.3,0.2,0.4,0.6");  // 003, ternary metal-oxide prototype (ICSD #36310)
        vparameters.push_back("6.5502,0.833715,0.76028213,136.1848,0.1867,0.3299,0.253,0.397,0.6382,0.226,0.897,0.6523");  // 004, ternary metal-oxide prototype (ICSD #182751)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mP40_14_2e_6e_2e"){
        vparameters.push_back("9.8340572172,0.944454808597,0.553242905653,103.024236752,0.2509,0.3482,0.0296,0.2556,0.9865,0.037,0.8687,0.659,0.6887,0.1227,0.5,0.7987,0.1052,0.724,0.4962,0.3748,0.159,0.2478,0.6341,0.016,0.2471,0.6018,0.303,0.1428,0.0441,0.6597,0.7519,0.5525,0.1628,0.3184");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_mC24_15_2e_af_e"){
        vparameters.push_back("6.0875108932,1.66963707964,1.11173475337,123.707427092,0.1695,0.483,0.8419,0.0299,0.6684,0.5221");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mC40_15_2e_3f_f"){
        vparameters.push_back("10.7421147996,0.897270659587,0.978465510971,147.31111283,0.45571,0.83777,0.78937,0.64875,0.16687,0.71802,0.50232,0.32358,0.43573,0.74259,0.03818,0.719,0.66021,0.26127");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_mC48_15_aef_3f_2e"){
        vparameters.push_back("5.8157729384,1.73589080063,2.00719678852,109.149281309,0.8352,0.16749,0.49971,0.9884,0.6705,0.0001,0.75396,0.48644,0.8628,0.21836,0.66556,0.86281,0.25504,0.34412,0.86546"); //DX20210427 - equivalent to Na2O3Pr (http://aflow.org/prototype-encyclopedia/A2B3C_mC48_15_aef_3f_2e.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC7_mC48_15_2f_e_e3f"){
        vparameters.push_back("13.2432547747,0.690228819766,0.966427058766,155.232510527,0.9409,0.2269,0.8889,0.6633,0.303,0.3712,0.8094,0.241,0.2004,0.6988,0.5659,0.2806,0.9947,0.1491,0.6373,0.8064,0.5797");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mC60_15_cf_e4f_ef"){
        vparameters.push_back("12.5843824384,0.579976673828,0.889919034503,68.755,0.6408,0.16,0.587,0.2597,0.0011,0.6275,0.0955,0.4035,0.729,0.1009,0.1137,0.6109,0.8071,0.2522,0.5497,0.5564,0.1061,0.6278,0.047,0.2555");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_33_a_a_2a"){
        vparameters.push_back("5.9547396971,1.2544327646,0.947018118432,0.93178,0.87342,-0.0,0.5759,0.8771,0.503,0.9576,0.9125,0.656,0.6123,0.839,0.9169");
        vparameters.push_back("5.301,1.0411243,1.3584229,0,0.061,0.123,0.006,0.59,0.874,0.064,0.395,0.173,0.172,0.94,0.898");  // 002, ternary metal-oxide prototype (ICSD #4416)
        vparameters.push_back("5.402,1.1795631,0.92687893,0.9179,0.8737,0,0.5793,0.8733,0.4936,0.5934,0.8612,0.8927,0.9303,0.8879,0.3708");  // 003, ternary metal-oxide prototype (ICSD #18152) //DX20210428 - equivalent to GaLiO2 (http://aflow.org/prototype-encyclopedia/ABC2_oP16_33_a_a_2a.html, part 3)
        vparameters.push_back("5.672,1.25811,0.94799013,0.925,0.872,0,0.6,0.865,0.5,0.95,0.9,0.33,0.62,0.875,0.9");  // 004, ternary metal-oxide prototype (ICSD #27117)
        vparameters.push_back("5.498,1.3106584,0.96362314,0.9367,0.874,0,0.4341,0.1231,0.0111,0.6153,0.8351,0.08,0.0364,0.0864,0.1584");  // 005, ternary metal-oxide prototype (ICSD #36652)
        vparameters.push_back("5.2728,1.2851426,1.0170308,0.9402,0.6427,0.9931,0.9294,0.1239,0.0013,0.9658,0.0679,0.6881,0.8883,0.6681,0.5654");  // 006, ternary metal-oxide prototype (ICSD #160643)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC24_36_b_ab_a"){
        vparameters.push_back("10.2764819278,0.576935497153,0.499471075467,0.108,0.155,0.167,0.504,0.33,0.83,0.0,0.353,0.81,0.594");
        vparameters.push_back("11.8784,0.57718211,0.4611648,0.077,0.895,0.166,0.563,0.666,0.839,0,0.63,0.786,0.5");
        vparameters.push_back("11.47,0.67323452,0.5603313,0.65,0.182,0.9111,0.25,0.1656,0.3571,0.25,0.3805,0.5,0");  // 003, ternary metal-oxide prototype (ICSD #1181)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_oP32_58_eg_3gh_g"){
        vparameters.push_back("7.7975902487,1.01536533044,0.715127549044,0.2419,0.8702,0.3614,0.922,0.1368,0.574,0.6392,0.8969,0.6001,0.754,0.7495,0.7697,0.8676,0.2387"); //DX20210427 - equivalent to Al2O5Si (SO_{2}, andalusite, http://aflow.org/prototype-encyclopedia/A2B5C_oP32_58_eg_3gh_g.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_oC28_63_c_ac_fg"){
        vparameters.push_back("6.2543235982,1.56175972927,1.21827411167,0.3472,0.7,0.75,0.5694,0.7778,0.45"); //DX20210428 - equivalent to CrNa2O4 (H1_{8}, http://aflow.org/prototype-encyclopedia/AB2C4_oC28_63_c_bc_fg.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_oC32_63_c_c2f_f"){
        vparameters.push_back("3.8117041359,2.59759679576,2.66755674233,0.1943,0.7775,0.9542,0.886,0.687,0.9362,0.8659,0.4343"); //DX20210427 - equivalent to Fe2O5Ti (pseudobrookite, http://aflow.org/prototype-encyclopedia/A2B5C_oC32_63_f_c2f_c.html, part 3)
        vparameters.push_back("3.6913,2.7012164,2.9846666,0.8878,0.3044,0.9576,0.6273,0.7628,0.4213,0.7983,0.5966");  // 002, ternary metal-oxide prototype (ICSD #50979)
        vparameters.push_back("3.6047,2.7507698,3.1203429,0.1008,0.7192,0.0471,0.6272,0.2449,0.4255,0.2057,0.6013");  // 003, ternary metal-oxide prototype (ICSD #50980)
        vparameters.push_back("3.415,3.55666471449,4.42079062958,0.69692,0.31968,0.58282,0.13527,0.1485,-0.05026,0.22082,0.10879");  // 004, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tI24_88_a_b_f"){
        vparameters.push_back("5.7543,2.2969084,0.8831,0.5138,0.2973");
        vparameters.push_back("5.7062,2.20627,0.2638,0.1105,0.4568"); //DX20210428 - equivalent to CaO4W (scheelite, H0_{4}, http://aflow.org/prototype-encyclopedia/AB4C_tI24_88_b_f_a.html, part 3)
        vparameters.push_back("5.213,2.1918281,0.9,0.5,0.325");  // 003, ternary metal-oxide prototype (ICSD #77334)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2_tP28_91_2d_b_ac"){
        vparameters.push_back("6.0709632596,1.40149875103,0.233,0.246,0.246,0.975,0.73,0.25,0.513,0.264,0.231");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_hP56_173_2b2c_ac_b5c"){
        vparameters.push_back("11.1609261164,0.841147428044,0.25,0.95,0.56,0.74,0.163,0.341,0.056,0.152,0.325,0.45,0.507,0.009,0.253,0.181,0.0,0.978,0.69,0.005,0.058,0.498,0.181,0.004,0.175,0.498,0.993,0.1161,0.32,0.255");
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_cF56_227_b_c_e"){
        vparameters.push_back("9.8805440428,0.2617");
      }
      // -------------------------------------------------------------------------
      // nitride prototypes (from R. Friedrich)
      // -------------------------------------------------------------------------
      // binaries
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_cI80_206_ad_e"){
        vparameters.push_back("9.8892799651,0.7716,0.1259,0.6002,0.6475");
      }
      // -------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_mC56_9_6a_6a_2a"){ //DX20210604 - added nitride for R. Friedrich (N3Na3W1, ICSD #75364) 
        vparameters.push_back("13.81,0.795293,0.887176,152.415,0.602,0.856,0.1,0.691,0.123,0.16,0.999,0.86,0.253,0.505,0.622,0.213,0.78,0.469,0.177,0.9,0.457,0.546,0.4729,0.5309,0.98,0.8144,0.733,0.436,0.81,0.733,0.185,0.708,0.002,0.971,0.4876,0.724,0.371,0.732,0.763,0.865,0.5,0.01238,0,0.681555,0.499645,0.4312");
      }
      if(anrl_label=="ABC_oP12_62_c_c_c"){
        vparameters.push_back("9.95851159,0.433951127369,0.653641836852,0.15617,0.97349,0.04761,0.42551,0.40749,0.75982");
        vparameters.push_back("7.1586,0.48988629,0.70044422,0.0948,0.4719,0.1442,0.9797,0.3792,0.731");  // 002, ternary metal-nitride prototype (ICSD #93259)
        vparameters.push_back("7.86,0.59287532,1.1119593,0.482,0.317,0.352,0.939,0.732,0.886");  // 003, metallic ternary prototype (ICSD #42757) //DX20210428 - equivalent to CuMnP (http://aflow.org/prototype-encyclopedia/ABC_oP12_62_c_c_c.html, part 3)
        vparameters.push_back("7.146,0.62244612,1.0720683,0.8048,0.9093,0.6959,0.589,0.012,0.2975");  // 004, metallic ternary prototype (ICSD #54301)
        vparameters.push_back("7.199,0.62911515,1.0247257,0.859,0.435,0.017,0.7973,0.204,0.371");  // 005, metallic ternary prototype (ICSD #54415)
        vparameters.push_back("7.77,0.61132561,1.0862291,0.731,0.901,0.494,0.296,0.387,0.917");  // 006, metallic ternary prototype (ICSD #58762)
        vparameters.push_back("6.625,0.65403774,1.1375094,0.198,0.431,0.462,0.715,0.7839,0.408");  // 007, metallic ternary prototype (ICSD #105338)
        vparameters.push_back("7.156,0.6352711,1.0444382,0.839,0.439,0.974,0.794,0.17,0.426");  // 008, metallic ternary prototype (ICSD #106459)
        vparameters.push_back("7.671,0.60643984,0.98970147,0.478,0.809,0.688,0.452,0.815,0.07");  // 009, metallic ternary prototype (ICSD #108571)
        vparameters.push_back("7.2813,0.62620686,1.0674605,0.8549771,0.067753334,0.96462031,0.67830856,0.70964818,0.38654465");  // 010, metallic ternary prototype (ICSD #150172)
        vparameters.push_back("7.7672,0.61037439,0.94261767,0.9694,0.3161,0.176,0.9257,0.3232,0.6007");  // 011, metallic ternary prototype (ICSD #157921)
        vparameters.push_back("6.842,0.59295528,1.1843028,0.7486,0.6244,0.9579,0.3169,0.8643,0.938");  // 012, metallic ternary prototype (ICSD #159305)
        vparameters.push_back("6.723,0.62977837,1.1518667,0.012,0.304,0.677,0.551,0.794,0.903");  // 013, metallic ternary prototype (ICSD #630573)
        vparameters.push_back("6.281,0.6938704,1.2023563,0.735,0.619,0.858,0.938,0.978,0.315");  // 014, metallic ternary prototype (ICSD #635080)
        vparameters.push_back("6.568,0.53562728,1.3532278,0.8957,0.5022,0.6879,0.0478,0.0038,0.8203");  // 015, ternary metal-oxide prototype (ICSD #65170)
        vparameters.push_back("8.75,0.42742857,0.74628571,0.19,0.739,0.374,0.529,0.094,0.858");  // 016, ternary metal-carbo-nitride prototype (ICSD #77321)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_tI10_139_e_e_a"){
        vparameters.push_back("4.0191421247,3.53238454024,0.664,0.8545");
        vparameters.push_back("3.8568,3.3538166,0.1449,0.3409");  // 002, ternary metal-nitride prototype (ICSD #80376)
        vparameters.push_back("4.152,3.1442678,0.6557,0.8589");  // 003, ternary metal-nitride prototype (ICSD #80377)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_hP6_191_f_a_d"){
        vparameters.push_back("5.3369,0.57059342");
      }

#if !(USE_HARDCODED_PROTOTYPES) //DX20210114 - these new prototypes are not hard-coded
      // -------------------------------------------------------------------------
      // metal-nitride prototypes (from DX)
      // ---------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_6_3ab_ab"){
        vparameters.push_back("4.86,0.58148148,0.84053498,99.49,0.192,0.064,0.825,0.654,0.62,0.359,0.201,0.425,0.353,0.805,0.778,0.028");  // 001, binary metal-nitride prototype (ICSD #187448)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC12_8_4a_2a"){
        vparameters.push_back("4.964,0.57836422,1.9284851,107.49,0.252,0.798,0.162,0.643,0.008,0.376,0.407,0.065,0.999,0.936,0.749,0.505");  // 001, binary metal-nitride prototype (ICSD #187447)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mC16_8_6a_b"){
        vparameters.push_back("9.0123,0.83308367,0.70781044,134.9427,0.86936937,0.5,0.12962963,0.5,0.63036963,0.76073926,0.36862937,0.23926074,0,0.4995005,0.5004995,0.5004995,0.5,0.75,0");  // 001, binary metal-nitride prototype (ICSD #155169)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_11_2e_e"){
        vparameters.push_back("5.1164,0.56651161,0.82880541,113.029,0.29511678,0.65490933,0.17069835,0.88921147,0.71680688,0.70665402");  // 001, binary metal-nitride prototype (ICSD #187444)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_mP14_11_e_6e"){
        vparameters.push_back("9.63,0.45794393,0.56386293,99.34,0.7819,0.8285,0.093,0.86,0.108,0.646,0.112,0.437,0.416,0.221,0.54,0.25,0.667,0.298");  // 001, binary metal-nitride prototype (ICSD #14244)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mC8_12_di_a"){
        vparameters.push_back("6.2233,0.58666624,0.85740684,71.556,0.0939,0.2687");  // 001, binary metal-nitride prototype (ICSD #29370)
        vparameters.push_back("6.3105,0.57963711,0.86886934,68.788,0.0836,0.2705");  // 002, binary metal-nitride prototype (ICSD #29376)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC12_12_2i_i"){
        vparameters.push_back("6.817,0.41587208,0.84264339,83.9924,0.373,0.864,0.902,0.47,0.675,0.769");  // 001, binary metal-nitride prototype (ICSD #187441)
        vparameters.push_back("5.2169,0.5773352,1.4982653,116.7975,0.33306,0,0.6369,0.45535,0.85,0.775");  // 002, binary metal-boride prototype (ICSD #418398)
        vparameters.push_back("7.2076,0.53114768,1.200247,125.55,0.6,0.021,0.468,0.42,0.0575,0.2544");  // 003, binary metal-carbide prototype (ICSD #54185)
        vparameters.push_back("7.2286,0.5329386,1.1971198,125.5426,0.626,0.065,0.522,0.447,0.04,0.2486");  // 004, binary metal-carbide prototype (ICSD #54188) //DX20210428 - equivalent to CaC2-III (http://aflow.org/prototype-encyclopedia/A2B_mC12_12_2i_i.html, part 3)
        vparameters.push_back("6.6749,0.43431362,1.0338582,120.885,0.399,0.112,0.857,0.469,0.154,0.2");  // 005, binary metal-carbide prototype (ICSD #181488)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC16_12_2i_2i"){
        vparameters.push_back("13.46,0.2830312,0.49946508,85.287,0.7559,0.7486,0.9779,0.0833,0.846,0.3926,0.6441,0.092");  // 001, binary metal-nitride prototype (ICSD #411555)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_mC30_12_a4i_3i"){
        vparameters.push_back("15.0498,0.25136547,0.62120427,79.46,0.62567,0.72225,0.328,0.96482,0.28973,0.62051,0.96826,0.66055,0.14105,0.51172,0.18256,0.86977,0.46936,0.81225");  // 001, binary metal-nitride prototype (ICSD #162794)
        vparameters.push_back("14.095,0.25902802,0.63213906,81.748,0.328,0.965,0.299,0.629,0.631,0.715,0.974,0.658,0.141,0.508,0.18,0.868,0.465,0.807");  // 002, binary metal-nitride prototype (ICSD #169726)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_13_g_e"){
        vparameters.push_back("5.053,0.58321789,0.98495943,129.9799,0.809,0.679,0.698,0.37");  // 001, binary metal-nitride prototype (ICSD #187451)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mP12_14_e_2e"){
        vparameters.push_back("5.6053,0.86029294,0.85947229,126.0785,0.2336,0,0.0131,0.1841,0.4,0.8837,0.3261,0.5888,0.169");  // 001, binary metal-nitride prototype (ICSD #160623)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mP16_14_e_3e"){
        vparameters.push_back("5.9607,2.1122855,0.60051672,103.253,0.2956,0.03505,0.3855,0.949,0.0987,0.19,0.914,0.1852,0.297,0.865,0.2681,0.362");  // 001, binary metal-nitride prototype (ICSD #98661)
        vparameters.push_back("4.4,1.47727272727,1.41869318182,117.78562,0.144,0.042,0.146,0.36,0.14,0.166,0.098,-0.005,0.294,0.877,0.112,-0.019");  // 002, (part 3)
        vparameters.push_back("4.46,1.94618834081,1.0201793722,120.5,0.3698,0.4285,0.0119,0.0865,0.4057,0.7725,0.4789,0.1446,0.7686,0.5981,0.4237,-0.0903");  // 003, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC12_15_e_f"){
        vparameters.push_back("7.1712,0.61281236,1.2249554,127.232,0.5516,0.491,0.607,0.032");  // 001, binary metal-nitride prototype (ICSD #280681) //DX20210428 - equivalent to C2Th (C_{g}, http://aflow.org/prototype-encyclopedia/A2B_mC12_15_f_e.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC16_36_2a_2a"){
        vparameters.push_back("2.59,4.6100386,1.7722008,0.86,0,0.77,0.82,0.47,0.82,0.15,0.93");  // 001, binary metal-nitride prototype (ICSD #167514)
        vparameters.push_back("5.0045,1.1481067,2.2077131,0.3496,0.9686,0.0567,0.8012,0,0,0.41054,0.77174");  // 002, binary metal-oxide prototype (ICSD #424729)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_53_h_h"){
        vparameters.push_back("2.79,1.344086,1.5197133,0.35,0.7,0.86,0.89");  // 001, binary metal-nitride prototype (ICSD #167513)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_57_e_c"){
        vparameters.push_back("2.901,1.7087211,2.7083764,0.297,0.792,0.08,0.838");  // 001, binary metal-nitride prototype (ICSD #187450)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_60_d_c"){
        vparameters.push_back("4.425,1.2481356,1.0915254,0.875,0.25,0.125,0.0833");  // 001, binary metal-nitride prototype (ICSD #150889) //DX20210428 - equivalent to zeta-Fe2N (http://aflow.org/prototype-encyclopedia/A2B_oP12_60_d_c.Fe2N.html, part 3)
        vparameters.push_back("4.515,1.2174972,1.0939092,0.829,0.714,0.376,0.588");  // 002, binary metal-oxide prototype (ICSD #15328)
        vparameters.push_back("4.948,1.2027082,1.1109539,0.822,0.74,0.43,0.56");  // 003, binary metal-oxide prototype (ICSD #20362) //DX20210428 - equivalent to alpha-PbO2 (http://aflow.org/prototype-encyclopedia/A2B_oP12_60_d_c.html, part 3)
        vparameters.push_back("4.8094,1.1733896,0.95660581,0.89,0.75,0.64,0.125");  // 004, binary metal-oxide prototype (ICSD #24060)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_oC32_63_c2f_cf"){
        vparameters.push_back("3.893,2.6365271,2.6365271,0.2373,0.70243,0.047,0.8812,0.3097,0.9268,0.13429,0.4399");  // 001, binary metal-nitride prototype (ICSD #16253)
        vparameters.push_back("3.8862,2.6277083,2.6407287,0.23678,0.8029,0.04701,0.88051,0.30862,0.92622,0.13455,0.44094");  // 002, binary metal-nitride prototype (ICSD #66533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC12_65_gj_e"){
        vparameters.push_back("5.701,1.3513419,0.49587792,0.117,0.31");  // 001, binary metal-nitride prototype (ICSD #187452)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI4_71_b_c"){
        vparameters.push_back("2.8796,1.479268,1.0460133");  // 001, binary metal-nitride prototype (ICSD #53146)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oI12_71_2e_f"){
        vparameters.push_back("13.0613,0.28741396,0.215415,0.45057355,0.20799023,0.84176015");  // 001, binary metal-nitride prototype (ICSD #187449)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI16_71_abe_n"){
        vparameters.push_back("10.7912,0.41118689,0.2889484,0.751,0.12321,0.6466");  // 001, binary metal-nitride prototype (ICSD #423831)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oI16_72_b_cj"){
        vparameters.push_back("5.94,0.94107744,1.0185185,0.855,0.145");  // 001, binary metal-nitride prototype (ICSD #27135)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B4_tI18_87_ah_h"){
        vparameters.push_back("6.873,0.62534556,0.6,0.2,0.9,0.3");  // 001, binary metal-nitride prototype (ICSD #26251)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI32_88_d_cf"){
        vparameters.push_back("8.653,0.64648099,0.173,0.827,0.625");  // 001, binary metal-nitride prototype (ICSD #30633) //DX20210428 - equivalent to Copper (I) Azide (CuN3, http://aflow.org/prototype-encyclopedia/AB3_tI32_88_d_cf.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP2_123_d_a"){
        vparameters.push_back("4.129,0.7808186");  // 001, binary metal-nitride prototype (ICSD #168645)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP4_123_ag_d"){
        vparameters.push_back("4.497,0.82432733,0.304");  // 001, binary metal-nitride prototype (ICSD #16963)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP8_123_egh_ab"){
        vparameters.push_back("3.745,2,0.25,0.25");  // 001, binary metal-nitride prototype (ICSD #180237)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP6_127_g_b"){
        vparameters.push_back("4.408,0.60528584,0.61441863");  // 001, binary metal-nitride prototype (ICSD #187442)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_131_e_c"){
        vparameters.push_back("2.704,2.2525888");  // 001, binary metal-nitride prototype (ICSD #187709)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_139_a_b"){
        vparameters.push_back("3.009,1.3921569");  // 001, binary metal-nitride prototype (ICSD #106932)
        vparameters.push_back("4.261,0.99835719");  // 002, binary metal-oxide prototype (ICSD #174027) //DX20210428 - equivalent to "Martensite Type" FeC_{x} )(L2_{0}, http://aflow.org/prototype-encyclopedia/AB_tI4_139_b_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_tI10_139_ae_e"){
        vparameters.push_back("4.2,2.8878571,0.667,0.8412");  // 001, binary metal-nitride prototype (ICSD #71638)
        vparameters.push_back("2.974,4.0773369,0.33341,0.1593");  // 002, binary metal-nitride prototype (ICSD #84202)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_tI14_139_ad_ce"){
        vparameters.push_back("4.382,1.9698768,0.2521");  // 001, binary metal-nitride prototype (ICSD #76389)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B_tI18_139_deh_a"){
        vparameters.push_back("5.72,1.0996503,0.293,0.757");  // 001, binary metal-nitride prototype (ICSD #41953) //DX20210428 - equivalent to Fe8N (D2_{g}, http://aflow.org/prototype-encyclopedia/A8B_tI18_139_deh_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_a_dh"){
        vparameters.push_back("6.106,1.1578775,0.635");  // 001, binary metal-nitride prototype (ICSD #24007)
        vparameters.push_back("6.5412,1.2368984,0.374");  // 002, binary metal-nitride prototype (ICSD #25008)
        vparameters.push_back("5.52,1.009058,0.84");  // 003, binary metal-nitride prototype (ICSD #183201)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_a_ch"){
        vparameters.push_back("6.1624,1.0959853,0.4359");  // 001, binary metal-nitride prototype (ICSD #187018)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI16_140_c_ah"){
        vparameters.push_back("6.1188,1.1606034,0.36406");  // 001, binary metal-nitride prototype (ICSD #1145)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI12_141_a_e"){
        vparameters.push_back("4.14,2.1268116,0.64");  // 001, binary metal-nitride prototype (ICSD #23403)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR4_160_3a_a"){
        vparameters.push_back("3.6458651,4.1802255,0.428,0.5,0.583,0");  // 001, binary metal-nitride prototype (ICSD #644523)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hR5_160_3a_2a"){
        vparameters.push_back("2.8499718,5.1435171,0.14085,0.68246,0.91112,0.79584,0.29163");  // 001, binary metal-nitride prototype (ICSD #185490)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_162_ek_ci"){
        vparameters.push_back("5.665,0.97440424,0.25,0.83333,0.494,0.75");  // 001, binary metal-nitride prototype (ICSD #43559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP5_164_ac_d"){
        vparameters.push_back("2.82,2.2163121,0.393,0.802");  // 001, binary metal-nitride prototype (ICSD #182700)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_164_ci_di"){
        vparameters.push_back("5.72,0.97902098,0.25,0,0.511,0.25,0.83333,0");  // 001, binary metal-nitride prototype (ICSD #60168)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hR3_166_a_c"){
        vparameters.push_back("3.8520969,5.3707781,0.2678");  // 001, binary metal-nitride prototype (ICSD #23530)
        vparameters.push_back("3.616977,4.9654486,0.7413");  // 002, binary metal-carbide prototype (ICSD #22283)
        vparameters.push_back("5.0166655,1.4456791,0.5867");  // 003, binary metal-carbide prototype (ICSD #186576)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_186_ac_bc"){
        vparameters.push_back("5.745,0.97859008,0.25,0,0.51,0.25,0.833,0");  // 001, binary metal-nitride prototype (ICSD #76280)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_187_a_h"){
        vparameters.push_back("3.2689,1.1635107,0.687679");  // 001, binary metal-nitride prototype (ICSD #290427)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP3_187_h_a"){
        vparameters.push_back("2.923,1.3373247,0.68184");  // 001, binary metal-nitride prototype (ICSD #290433)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_189_f_ad"){
        vparameters.push_back("5.196,0.56023865,0.3928");  // 001, binary metal-nitride prototype (ICSD #1396)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP3_191_e_a"){
        vparameters.push_back("2.829,1.7575115,0.384");  // 001, binary metal-nitride prototype (ICSD #260545)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_191_af_d"){
        vparameters.push_back("5.263,0.5166255");  // 001, binary metal-nitride prototype (ICSD #182351)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_191_bc_ad"){
        vparameters.push_back("3.98,1.2084925");  // 001, binary metal-nitride prototype (ICSD #161079)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_193_g_b"){
        vparameters.push_back("7.6418,0.9225706,0.2721");  // 001, binary metal-nitride prototype (ICSD #77730)
        vparameters.push_back("8.78,0.85649203,0.75");  // 002, binary metal-oxide prototype (ICSD #15695)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP4_194_a_c"){
        vparameters.push_back("3.31,1.7492447");  // 001, binary metal-nitride prototype (ICSD #185566)
        vparameters.push_back("3.358,1.2084574");  // 002, binary metal-boride prototype (ICSD #24363)
        vparameters.push_back("4.362,0.575424117377");  // 003, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_194_c_e"){
        vparameters.push_back("3.2647,2.3412565,0.5930479");  // 001, binary metal-nitride prototype (ICSD #290428)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_194_e_c"){
        vparameters.push_back("3.0543,2.549946,0.0905853");  // 001, binary metal-nitride prototype (ICSD #290431)
        vparameters.push_back("2.88,2.7777778,0.589");  // 002, binary metal-carbide prototype (ICSD #168280)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_194_bf_ac"){
        vparameters.push_back("2.8413,3.4114666,0.075");  // 001, binary metal-nitride prototype (ICSD #25656)
        vparameters.push_back("3.8911,3.3888361,0.5862");  // 002, binary metal-nitride prototype (ICSD #162797)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_194_cf_f"){
        vparameters.push_back("2.89,5.2892734,0.07,0.84");  // 001, binary metal-nitride prototype (ICSD #186207)
        vparameters.push_back("2.9051,4.4103473,0.4722,0.63985");  // 002, binary metal-boride prototype (ICSD #23715)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP16_194_bh_ag"){
        vparameters.push_back("5.72,0.96153846,0.52");  // 001, binary metal-nitride prototype (ICSD #106926)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_cI32_204_c_g"){
        vparameters.push_back("6.592,0.6933,0.0944");  // 001, binary metal-nitride prototype (ICSD #162105)
        vparameters.push_back("6.289,0.8855,0.3784");  // 002, binary metal-nitride prototype (ICSD #162106)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cP5_221_ac_b"){
        vparameters.push_back("3.868");  // binary metal-nitride prototype (ICSD #44369)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_cP5_221_a_bc"){
        vparameters.push_back("3.74");  // binary metal-nitride prototype (ICSD #76403) //DX20210428 - equivalent to Fe4N (http://aflow.org/prototype-encyclopedia/A4B_cP5_221_bc_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cP7_221_ac_d"){
        vparameters.push_back("4.122");  // binary metal-nitride prototype (ICSD #30370)
      }

      // -------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_aP8_1_a_2a_5a"){
        vparameters.push_back("6.57,0.6173516,0.56754947,90.52,90.53,91.16,0.788,0,0.759,0.791,0.497,0.745,0.291,0.505,0.238,0.958,0.555,0.254,0.124,0.447,0.732,0.285,0.998,0.127,0.621,0.465,0.246,0.459,0.536,0.74");  // 001, ternary metal-nitride prototype (ICSD #92316)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_mC16_5_a_c_a2bc"){
        vparameters.push_back("6.731,0.88307829,0.94829892,88.82,0.093,0.631,0.084,0.617,0.73,0.85,0.763,0.715,0.848,0.232");  // 001, ternary metal-nitride prototype (ICSD #92315)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C3_mP8_6_2ab_ab_a2b"){
        vparameters.push_back("6.32,0.58164557,0.60981013,90.31,0.972,0.462,0.283,0,0.297,0.5,0.638,0.566,0.819,0,0.814,0.502,0.473,0.395,0.143,0.617");  // 001, ternary metal-nitride prototype (ICSD #92312)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_mP8_6_3a2b_ab_b"){
        vparameters.push_back("6.6606,0.57952737,0.64258475,91.4279,0.363,0.521,0.67,0.005,0.026,0.447,0.688,0.505,0.221,0,0.535,0.61,0.207,0.498,0.872,0.659");  // 001, ternary metal-nitride prototype (ICSD #92314)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_mC20_12_2i_i_2i"){
        vparameters.push_back("10.93,0.45379689,0.86336688,141.6047,0.726,0.7657,0.0848,0.4066,0.5029,0.8751,0.7185,0.1958,0.4256,0.613");  // 001, ternary metal-nitride prototype (ICSD #72389)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C12_mC30_12_a_i_6i"){
        vparameters.push_back("14.272,0.26534473,0.62268778,87.17,0.347,0.6745,0.3979,0.0033,0.3277,0.0835,0.2593,0.1569,0,0.7378,0.068,0.6555,0.133,0.5702");  // 001, ternary metal-nitride prototype (ICSD #31297)
        vparameters.push_back("14.254,0.2655395,0.62319349,87.1,0.847,0.674,0.398,0.003,0.328,0.084,0.259,0.157,0,0.738,0.068,0.656,0.133,0.57");  // 002, ternary metal-nitride prototype (ICSD #659621)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mP12_14_e_e_e"){
        vparameters.push_back("4.54,1.0484581,1.2790749,124.9,0.91998,0.39152,0.608,0.67587,0.54156,0.91072,0.84932,0.19521,0.8007");  // 001, ternary metal-nitride prototype (ICSD #402341)
        vparameters.push_back("5.7612,0.986617371381,1.00107616469,111.721,0.14746,0.13055,0.86937,0.28353,-0.00643,0.29429,0.6551,0.1311,0.3211");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_2c_c_c"){
        vparameters.push_back("3.532,5.7202718,1.4068516,0.78124,0.62006,0.96137,0.19906");  // 001, ternary metal-nitride prototype (ICSD #96228)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC20_63_c_a_cf"){
        vparameters.push_back("3.369,3.4128822,2.6662214,0.2526,0.9554,0.3719,0.958");  // 001, ternary metal-nitride prototype (ICSD #29521)
        vparameters.push_back("3.1345,3.150646,2.3273887,0.74981,0.57459,0.12692,0.94936");  // 002, ternary metal-oxide prototype (ICSD #159026)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3_oC28_63_cg_c_cg"){
        vparameters.push_back("8.503,1.2094555,0.59179113,0.8921,0.3055,0.1243,0.2843,0.1174,0.6918,0.8714");  // 001, ternary metal-nitride prototype (ICSD #40205)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC8_65_i_c_b"){
        vparameters.push_back("3.698,1.7212006,1.2439156,0.166");  // 001, ternary metal-nitride prototype (ICSD #92307)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C3_oC24_65_a2ij_i_cj"){
        vparameters.push_back("4.9,5.0531224,0.99591837,0.1982,0.40004,0.1027,0.0998,0.29317");  // 001, ternary metal-nitride prototype (ICSD #98178)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_oI12_71_e_a_de"){
        vparameters.push_back("12.4601,0.30095264,0.2764424,0.35332,0.1638");  // 001, ternary metal-nitride prototype (ICSD #50579)
        vparameters.push_back("12.9655,0.31627781,0.30133045,0.64,0.842");  // 002, ternary metal-oxide prototype (ICSD #68217)
        vparameters.push_back("18.3959,0.1845683,0.1801434,0.81646899,0.68969595");  // 003, ternary metal-oxide prototype (ICSD #93651)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_99_ab_ab_ab"){
        vparameters.push_back("3.0285,2.430246,0.76,0.05,0.35,0.58,0.3,0");  // 001, ternary metal-nitride prototype (ICSD #23779)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C3_tP8_115_abc_g_de"){
        vparameters.push_back("3.8929,2.021886,0.23552817,0.098179023");  // 001, ternary metal-nitride prototype (ICSD #92311)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP3_123_c_a_b"){
        vparameters.push_back("2.83,1.3120141");  // 001, ternary metal-nitride prototype (ICSD #53505)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_tP8_123_ai_bc_d"){
        vparameters.push_back("3.965,1.3881463,0.75");  // 001, ternary metal-nitride prototype (ICSD #92313)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_129_c_c_c"){
        vparameters.push_back("3.5701,2.1170836,0.65036,0.08393,0.3303");  // 001, ternary metal-nitride prototype (ICSD #2027)
        vparameters.push_back("3.081,2.5316456,0.9,0.617,0.335");  // 002, ternary metal-nitride prototype (ICSD #100437)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_129_c_c_bc"){
        vparameters.push_back("4.1279,2.0304755,0.8479,0.4142,0.168");  // 001, ternary metal-nitride prototype (ICSD #50994)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_129_a_c_bc"){
        vparameters.push_back("3.895,1.5697047,0.221,0.789");  // 001, ternary metal-nitride prototype (ICSD #92309)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_131_e_c_b"){
        vparameters.push_back("3.5809,1.9574967");  // 001, ternary metal-nitride prototype (ICSD #69044)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tP6_131_b_c_e"){
        vparameters.push_back("3.924,1.8055556");  // 001, ternary metal-nitride prototype (ICSD #87414)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI8_139_d_b_a"){
        vparameters.push_back("3.653,1.460991");  // 001, ternary metal-nitride prototype (ICSD #92305)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B6C2_tI22_139_ae_eg_e"){
        vparameters.push_back("3.905,5.1766965,0.1805,0.6949,0.6023,0.5871");  // 001, ternary metal-nitride prototype (ICSD #98477)
        vparameters.push_back("4.0698,4.9520861,0.82019,0.3019,0.39895,0.9153");  // 002, ternary metal-nitride prototype (ICSD #411473)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_tI20_140_h_h_a"){
        vparameters.push_back("5.6213,1.3433192,0.872,0.333");  // 001, ternary metal-nitride prototype (ICSD #413356)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI20_140_h_a_h"){
        vparameters.push_back("5.5615,1.2370224,0.3794,0.8312");  // 001, ternary metal-nitride prototype (ICSD #413357)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tI20_140_a_h_h"){
        vparameters.push_back("5.6646,1.4818699,0.372,0.836");  // 001, ternary metal-nitride prototype (ICSD #415304)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR4_160_a_2a_a"){
        vparameters.push_back("2.8560976,5.4641029,0.8251,0.2629,0.4068,0");  // 001, ternary metal-nitride prototype (ICSD #84639)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_hP5_164_d_d_a"){
        vparameters.push_back("3.282,1.6636197,0.39,0.77");  // 001, ternary metal-nitride prototype (ICSD #16231)
        vparameters.push_back("3.622,1.7556599,0.6314,0.2797");  // 002, ternary metal-nitride prototype (ICSD #410826) and Ce2O2S (http://aflow.org/prototype-encyclopedia/A2B2C_hP5_164_d_d_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_hP5_164_a_d_d"){
        vparameters.push_back("3.557,1.5451223,0.65,0.27");  // 001, ternary metal-nitride prototype (ICSD #34003)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR4_166_a_c_b"){
        vparameters.push_back("3.1359829,5.5605955,0.881");  // 001, ternary metal-nitride prototype (ICSD #71136)
        vparameters.push_back("4.3180071,3.9173242,0.189");  // 002, metallic ternary prototype (ICSD #615828)
        vparameters.push_back("3.8272968,4.9635547,0.107");  // 003, ternary metal-oxide prototype (ICSD #165068)
        vparameters.push_back("3.9732098,3.7823428,0.082");  // 004, ternary metal-carbo-nitride prototype (ICSD #59860)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3_hP14_176_h_c_h"){
        vparameters.push_back("8.0141,0.69977914,0.35834333,0.086966667,0.12945233,0.68896467");  // 001, ternary metal-nitride prototype (ICSD #36502)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_b_b_a"){
        vparameters.push_back("3.3306,3.2702516,0,0.598,0.2466");  // 001, ternary metal-nitride prototype (ICSD #172471)
        vparameters.push_back("4.3972,1.3683253,0,0.777,0.245");  // 002, metallic ternary prototype (ICSD #54657)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_186_a_ab_b"){
        vparameters.push_back("2.9213,3.7507274,0,0.2,0.4,0.8");  // 001, ternary metal-nitride prototype (ICSD #80029)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP3_187_a_c_d"){
        vparameters.push_back("3.64,0.93956044");  // 001, ternary metal-nitride prototype (ICSD #247028)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C3_hP11_189_dg_g_f"){
        vparameters.push_back("6.475,0.54903475,0.3557,0.7281,0.3537");  // 001, ternary metal-nitride prototype (ICSD #411152)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hP4_191_c_a_b"){
        vparameters.push_back("3.65,1.260274");  // 001, ternary metal-nitride prototype (ICSD #92308)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP4_191_b_a_c"){
        vparameters.push_back("4,1.05");  // 001, ternary metal-nitride prototype (ICSD #92310)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP8_194_c_a_f"){
        vparameters.push_back("3.6506,3.4680053,0.595");  // 001, ternary metal-nitride prototype (ICSD #74791)
        vparameters.push_back("3.039,4.0786443,0.0833");  // 002, ternary metal-oxide prototype (ICSD #2786) //DX20210428 - equivalent to hexagonal delafossite (AlCuO2, http://aflow.org/prototype-encyclopedia/ABC2_hP8_194_a_c_f.html, part 3)
        vparameters.push_back("3.672,3.7363834,0.93");  // 003, ternary metal-oxide prototype (ICSD #27336)
        vparameters.push_back("3.035,3.7723229,0.0892");  // 004, ternary metal-oxide prototype (ICSD #66546)
        vparameters.push_back("3.04,4.9694079,0.5794");  // 005, ternary metal-nitride prototype (ICSD #165101)
        vparameters.push_back("3.151,4.427166,0.092");  // 006, ternary metal-nitride prototype (ICSD #181247)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_a_f_b"){
        vparameters.push_back("2.8724,3.8201504,0.128");  // 001, ternary metal-nitride prototype (ICSD #75971)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP8_194_a_c_f"){
        vparameters.push_back("2.9108,3.6245362,0.6384");  // 001, ternary metal-nitride prototype (ICSD #185913)
        vparameters.push_back("2.958,3.9148073,0.6345");  // 002, ternary metal-oxide prototype (ICSD #29282)
        vparameters.push_back("2.8869,4.2823097,0.075");  // 003, ternary metal-oxide prototype (ICSD #95663)
        vparameters.push_back("3.209,4.3035213,0.58");  // 004, ternary metal-nitride prototype (ICSD #42926)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_hP10_194_h_a_c"){
        vparameters.push_back("8.4414,0.82707845,0.8561");  // 001, ternary metal-nitride prototype (ICSD #67497)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_hP16_194_c_ae_2f"){
        vparameters.push_back("2.795,8.5144902,0.6078,0.0543,0.6591");  // 001, ternary metal-nitride prototype (ICSD #181351)
        vparameters.push_back("3.091,7.6700097,0.11125,0.05524,0.66016");  // 002, ternary metal-nitride prototype (ICSD #157843)
        vparameters.push_back("3.0935,8.0697915,0.10979,0.05448,0.65873");  // 003, ternary metal-nitride prototype (ICSD #159456)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cF32_227_b_c_a"){
        vparameters.push_back("8.74");  // ternary metal-nitride prototype (ICSD #72546)
      }

      // -------------------------------------------------------------------------
      // metal prototypes (from DX) //DX20201028
      // ---------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC16_8_b_2a_2ab"){
        vparameters.push_back("7.171,1.2100126,0.7774369,71.74,0.12,0,0.641,0.978,0.961,0.538,0.551,0.528,0.281,0.224,0.42,0.417,0.756,0.876");  // 001, metallic ternary prototype (ICSD #58830)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C3_mP9_10_m_2n_am"){
        vparameters.push_back("8.225,0.45483283,0.91671733,109.65,0.89515,0.35185,0.4983,0.7077,0.1823,0.39514,0.2748,0.90637");  // 001, metallic ternary prototype (ICSD #417035)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C4_mP11_10_m_a2m_2n"){
        vparameters.push_back("9.046,0.36325448,0.82367897,94.53,0.6346,0.5766,0.8578,0.343,0.6676,0.9576,0.1541,0.3154,0.5619,0.295");  // 001, metallic ternary prototype (ICSD #55578)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_mP10_11_e_2e_2e"){
        vparameters.push_back("10.783,0.45015302,0.45618103,91.3,0.7587,0.2575,0.4999,0.7506,0.1184,0.2145,0,0.6962,0.3728,0.2462");  // 001, metallic ternary prototype (ICSD #416299)
        vparameters.push_back("7.9823,0.36573419,0.70548589,108.493,0.2366,0.9736,0.1236,0.3867,0.5582,0.3219,0.6818,0.7651,0.8254,0.2724");  // 002, ternary metal-oxide prototype (ICSD #422751)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_mP10_11_2e_e_2e"){
        vparameters.push_back("11.3765,0.40427196,0.40718147,90.542,0.13236,0.2261,0.49959,0.7565,0.75464,0.2577,0.3606,0.2566,0.001,0.7307");  // 001, metallic ternary prototype (ICSD #424108)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mP14_11_2e_4e_e"){
        vparameters.push_back("9.591,0.47361068,0.79799812,107.838,0.05956,0.27432,0.57956,0.27679,0.8773,0.90331,0.57123,0.89518,0.36827,0.46211,0.87129,0.49589,0.21178,0.8337");  // 001, metallic ternary prototype (ICSD #261042)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_mC10_12_i_a_i"){
        vparameters.push_back("11.1,0.38801802,0.95268468,154.4287,0.5942,0.8555,0.9238,0.5702");  // 001, metallic ternary prototype (ICSD #182050)
        vparameters.push_back("5.0674,0.57735722,1.718475,112.8266,0.43068,0.64178,0.2491,0.8712");  // 002, ternary metal-oxide prototype (ICSD #172559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mC12_12_i_i_i"){
        vparameters.push_back("11.098,0.3973689,0.98442963,154.8009,0.2238,0.8102,0.8046,0.675,0.3393,0.1465");  // 001, metallic ternary prototype (ICSD #20632)
        vparameters.push_back("12.87,0.28135198,0.98036519,151.4655,0.218,0.346,0.492,0.24,0.787,0.149");  // 002, ternary metal-oxide prototype (ICSD #1570)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C2_mC14_12_ai_i_i"){
        vparameters.push_back("11.019,0.38767583,0.94896089,136.7589,0.30129,0.44343,0.42164,0.79411,0.8673,0.22347");  // 001, metallic ternary prototype (ICSD #410985)
        vparameters.push_back("10.47,0.32139446,0.99136581,149.1835,0.9531,0.237,0.27014,0.1147,0.70952,0.31653");  // 002, ternary metal-nitride prototype (ICSD #62083)
        vparameters.push_back("11.616,0.2858385,0.93737087,151.0914,0.4779,0.7278,0.2203,0.8301,0.82,0.6689");  // 003, ternary metal-nitride prototype (ICSD #88511)
        vparameters.push_back("12.0319,0.2773793,0.91653854,151.5643,0.4875,0.7375,0.79392,0.18811,0.17052,0.32796");  // 004, ternary metal-nitride prototype (ICSD #417827)
        vparameters.push_back("11.563,0.28629248,0.94049122,151.087,0.961,0.236,0.714,0.324,0.272,0.12");  // 005, ternary metal-nitride prototype (ICSD #617669)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_mC18_12_3i_i_a"){
        vparameters.push_back("10.077,0.54314776,0.92448149,81.93,0.5297,0.2753,0.7074,0.9295,0.8264,0.5908,0.85526,0.24903");  // 001, metallic ternary prototype (ICSD #260159)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C4_mC22_12_i_a2i_2i"){
        vparameters.push_back("16.9369,0.27129522,0.97586926,153.6114,0.38824,0.6115,0.04121,0.23355,0.40758,0.34101,0.71885,0.06875,0.2201,0.643");  // 001, metallic ternary prototype (ICSD #419134)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C8_mC30_12_ai_hi_2ij"){
        vparameters.push_back("18.549,0.36039679,0.36616529,81.66,0.213,0.67755,0.6913,0.791,0.13,0.8708,0.4541,0.51517,0.202,0.14945,0.26752,0.19789");  // 001, metallic ternary prototype (ICSD #240016)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C7_mC30_12_3i_i_a3i"){
        vparameters.push_back("14.257,0.32012345,0.86806481,86.021,0.63732,0.49687,0.46213,0.20062,0.19481,0.20687,0.83052,0.11118,0.81799,0.32444,0.01039,0.38194,0.65953,0.01844");  // 001, metallic ternary prototype (ICSD #171243)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_mC20_15_e_f_f"){
        vparameters.push_back("10.487,0.5754744,0.79202823,76.32,0.1345,0.63313,0.36433,0.99699,0.85144,0.3604,0.64786");  // 001, metallic ternary prototype (ICSD #391432)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mC28_15_f_2f_e"){
        vparameters.push_back("19.0003,0.35246812,0.35005237,86.065,0.7786,0.2981,0.0987,0.1999,0.42255,0.20721,0.77201,0.34791,0.51916,0.01422");  // 001, metallic ternary prototype (ICSD #261117)
        vparameters.push_back("9.7551,0.60957858,0.51185534,89.437,0.1702,0.674,0.138,0.195,0.64,0.876,0.379,0.889,0.08,0.563");  // 002, ternary metal-oxide prototype (ICSD #1044)
        vparameters.push_back("9.753,0.61047883,0.51204758,89.42,0.67163,0.201,0.147,0.246,0.6073,0.8998,0.4415,0.3737,0.358,0.1221");  // 003, ternary metal-oxide prototype (ICSD #14196)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_20_b_a_c"){
        vparameters.push_back("8.38,1.2601432,0.75536993,0.7698,0.504,0.351,0.673,0.513");  // 001, metallic ternary prototype (ICSD #1156)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_33_a_a_a"){
        vparameters.push_back("7.5356,1.0569563,0.62430331,0.51211,0.79877,0.75,0.29931,0.41498,0.7559,0.18658,0.0868,0.7533");  // 001, metallic ternary prototype (ICSD #106418)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_36_a_a_a"){
        vparameters.push_back("4.8891,1.720603,1.6021354,0.66671918,0.026102028,0.99999382,0.23452562,0.66654623,0.45737235");  // 001, metallic ternary prototype (ICSD #55819)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC32_36_b_2a_2ab"){
        vparameters.push_back("8.681,0.82363783,1.2135699,0.638,0.011,0.12,0,0.082,0.731,0.673,0.736,0.7763,0.622,0.29,0.756,0.619,0.563");  // 001, metallic ternary prototype (ICSD #56278)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC18_38_ad_be_bd"){
        vparameters.push_back("4.0307,3.2806709,1.8867442,0.74976771,0.40946388,0.0012150248,0.8747892,0.1251434,0.66729312,0.00026107832,0.70531985,0.29437221");  // 001, metallic ternary prototype (ICSD #59430)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10BC3_oP14_47_2q2rs_h_bt"){
        vparameters.push_back("4.261,1.0073926,3.6026754,0.3187,0.1603,0.0819,0.4156,0.4159,0.255");  // 001, metallic ternary prototype (ICSD #20663)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP8_51_f_e_be"){
        vparameters.push_back("5.254,0.78625809,1.3644842,0.2824,0.8886,0.656");  // 001, metallic ternary prototype (ICSD #106800)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP8_51_ae_f_f"){
        vparameters.push_back("5.09,0.79980354,1.3709234,0.315,0.118,0.683");  // 001, metallic ternary prototype (ICSD #623100)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_oP12_51_afj_e_e"){
        vparameters.push_back("7.7097,0.52203069,0.91056721,0.17610349,0.61101001,0.95660968,0.92790085,0.71070497");  // 001, metallic ternary prototype (ICSD #9986)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C_oP12_51_ak_ef_f"){
        vparameters.push_back("5.66,1.1909894,1.5980565,0.3741,0.1328,0.57444,0.20949,0.78478");  // 001, metallic ternary prototype (ICSD #249924)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP12_51_e_afj_e"){
        vparameters.push_back("8.676,0.48744813,0.85892116,0.5998,0.19219,0.93475,0.93556,0.68673");  // 001, metallic ternary prototype (ICSD #410891)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP12_57_c_d_d"){
        vparameters.push_back("5.75,1.3547826,0.9773913,0.104,0.6363,0.3994,0.8726,0.0244");  // 001, metallic ternary prototype (ICSD #370036)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_57_d_2d_c"){
        vparameters.push_back("6.804,1.96796,0.95517343,0.093,0.6601,0.8945,0.841,0.4885,0.611,0.1664");  // 001, metallic ternary prototype (ICSD #107616)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP8_59_a_e_b"){
        vparameters.push_back("4.216,1.2018501,1.0661765,0.333,0.333,0.5,0.833");  // 001, metallic ternary prototype (ICSD #628580)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_oP14_59_ae_ae_b"){
        vparameters.push_back("4.5526,1.7043008,1.9876554,0.04449,0.36673,0.2694,0.52383,0.58763,0.55436,0.88559");  // 001, metallic ternary prototype (ICSD #245679)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oP16_59_b_ae_ef"){
        vparameters.push_back("6.938,1.8080138,0.94018449,0.4367,0.68751,0.9941,0.1706,0.62693,0.33466,0.52489,0.95941");  // 001, metallic ternary prototype (ICSD #107444)
        vparameters.push_back("7.333,1.8847675,0.93113323,0.4473,0.6873,0.99594,0.1735,0.63365,0.33473,0.52812,0.962");  // 002, metallic ternary prototype (ICSD #107448)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_62_c_d_c"){
        vparameters.push_back("7.946,1.0187516,0.69808709,0.0332,0.89792,0.37781,0.11849,0.82793,0.9443,0.59466");  // 001, metallic ternary prototype (ICSD #425509)
        vparameters.push_back("7.38,0.95799458,0.75880759,0.379,0.652,0.033,0.35,0.823,0.947,0.088");  // 002, metallic ternary prototype (ICSD #602790)
        vparameters.push_back("5.307,1.2504051,0.83389862,0.39,0.55,0.036,0.352,0.814,0.937,0.172");  // 003, ternary metal-boride prototype (ICSD #603579)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_63_c_c_c"){
        vparameters.push_back("3.6684,4.8307164,1.086441,0.7978,0.95664,0.58754");  // 001, metallic ternary prototype (ICSD #166874)
        vparameters.push_back("3.212,4.3539851,0.96575342,0.696,0.535,0.9109");  // 002, ternary metal-boride prototype (ICSD #16777)
        vparameters.push_back("3.268,4.1156671,0.95318237,0.54,0.901,0.707");  // 003, ternary metal-boride prototype (ICSD #43005)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_63_c_c_a"){
        vparameters.push_back("4.2145,2.2335271,1.2607189,0.5972146,0.29158002");  // 001, metallic ternary prototype (ICSD #424600)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_f_c_c"){
        vparameters.push_back("4.08,2.4877451,1.7303922,0.213,0.939,0.348,0.956");  // 001, metallic ternary prototype (ICSD #57645)
        vparameters.push_back("4,2.3075,1.785,0.222,0.928,0.356,0.944");  // 002, metallic ternary prototype (ICSD #57693) //DX20210428 - equivalent to Al2CuMg (E1_{1}, http://aflow.org/prototype-encyclopedia/A2BC_oC16_63_f_c_c.html, part 3)
        vparameters.push_back("4.037,2.3661135,1.6019321,0.796,0.085,0.364,0.05");  // 003, metallic ternary prototype (ICSD #103875)
        vparameters.push_back("4.3642,2.589478,1.8801613,0.77214,0.05365,0.66459,0.94045");  // 004, metallic ternary prototype (ICSD #425493)
        vparameters.push_back("4.223,2.6381719,1.6663509,0.561,0.287,0.152,0.956");  // 005, metallic ternary prototype (ICSD #658140)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_a_c_g"){
        vparameters.push_back("7.7183,0.90569167,0.70617882,0.6485,0.819,0.3032");  // 001, metallic ternary prototype (ICSD #261785)
        vparameters.push_back("9.1362,0.8928548,0.62475646,0.3044,0.8289,0.7637");  // 002, metallic ternary prototype (ICSD #261786)
        vparameters.push_back("9.2835,0.8182582,0.63188453,0.3153,0.8185,0.6763");  // 003, metallic ternary prototype (ICSD #261788)
        vparameters.push_back("10.632,0.74473288,0.62046652,0.72594,0.8269,0.328");  // 004, metallic ternary prototype (ICSD #380341)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_c_2c"){
        vparameters.push_back("4.75,4.1263158,0.97894737,0.609,0.824,0.25,0.957");  // 001, metallic ternary prototype (ICSD #58647)
        vparameters.push_back("4.3923,3.6935091,0.98661294,0.6994,0.038,0.5648,0.2504");  // 002, metallic ternary prototype (ICSD #240098)
        vparameters.push_back("4.499,3.9319849,1,0.894,0.685,0.55,0.25");  // 003, metallic ternary prototype (ICSD #621687)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_63_c_f_c"){
        vparameters.push_back("4.2266,2.4376567,1.9779492,0.55971,0.27339,0.1543,0.9448");  // 001, metallic ternary prototype (ICSD #96152)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_g_c_c"){
        vparameters.push_back("7.097,1.5019022,0.6299845,0.6747,0.13117,0.7127,0.90293");  // 001, metallic ternary prototype (ICSD #99139)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC20_63_g_2c_a"){
        vparameters.push_back("7.0102,2.1704231,0.98316738,0.584,0.205,0.816,0.888");  // 001, metallic ternary prototype (ICSD #658701)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_oC24_63_acf_c_c"){
        vparameters.push_back("4.056,3.7302761,1.6346154,0.557,0.379,0.729,0.814,0.946");  // 001, metallic ternary prototype (ICSD #57760)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_oC24_63_c_cg_e"){
        vparameters.push_back("9.765,0.80747568,0.52892985,0.329,0.948,0.202,0.194,0.708");  // 001, metallic ternary prototype (ICSD #58898)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC3_oC28_63_cf_a_cf"){
        vparameters.push_back("4.098,2.4670571,3.1747194,0.3088,0.9983,0.5714,0.9114,0.289,0.8959");  // 001, metallic ternary prototype (ICSD #10044)
        vparameters.push_back("3.173,2.6542704,3.3810274,0.7047,0,0.5929,0.1136,0.2899,0.1064");  // 002, ternary metal-boride prototype (ICSD #25753)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC32_63_ac_f_2cf"){
        vparameters.push_back("4.639,3.6861393,2.0773874,0.7976,0.5921,0.4277,0.3515,0.4931,0.1961,0.4111");  // 001, metallic ternary prototype (ICSD #410732)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC10_65_j_i_a"){
        vparameters.push_back("4.0172,3.4774221,0.9045106,0.1986,0.3599");  // 001, metallic ternary prototype (ICSD #54612)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC12_65_i_cf_a"){
        vparameters.push_back("5.3308,1.7232498,0.69906956,0.66645548");  // 001, metallic ternary prototype (ICSD #107864)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C2_oC14_65_aj_j_i"){
        vparameters.push_back("3.9788,5.3289937,0.91944807,0.17641,0.4198,0.28698");  // 001, metallic ternary prototype (ICSD #240761)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_65_aci_j_i"){
        vparameters.push_back("4.192,4.1898855,0.98592557,0.288,0.416,0.136");  // 001, metallic ternary prototype (ICSD #103850)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_oC18_65_j_a_2ij"){
        vparameters.push_back("4.313,5.2124275,1.0159286,0.2046,0.4311,0.3222,0.0697");  // 001, metallic ternary prototype (ICSD #165783)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C7_oC24_65_aj_i_c2ij"){
        vparameters.push_back("4.565,5.9824754,1,0.875,0.593,0.7813,0.6847,0.9107");  // 001, metallic ternary prototype (ICSD #102239)
        vparameters.push_back("4.579,5.7896921,1,0.8718,0.592,0.7812,0.6845,0.9087");  // 002, metallic ternary prototype (ICSD #160913)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oI10_71_f_h_a"){
        vparameters.push_back("8.5967,0.66765154,0.51107983,0.7969,0.2835");  // 001, metallic ternary prototype (ICSD #55495)
        vparameters.push_back("8.994,0.60973983,0.47520569,0.1961,0.2537");  // 002, metallic ternary prototype (ICSD #58885)
        vparameters.push_back("8.523,0.71441981,0.36595096,0.184,0.225");  // 003, ternary metal-oxide prototype (ICSD #6158)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3_oI14_71_e_f_af"){
        vparameters.push_back("19.795,0.24976004,0.22556201,0.4351,0.8584,0.6832");  // 001, metallic ternary prototype (ICSD #9564)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_oI16_71_an_e_d"){
        vparameters.push_back("9.562,0.73457436,0.41612633,0.758,0.354,0.31");  // 001, metallic ternary prototype (ICSD #58046)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C4_oI22_71_n_af_eg"){
        vparameters.push_back("16.5484,0.51916802,0.27410505,0.80940162,0.81788308,0.65938423,0.61795128,0.69904145");  // 001, metallic ternary prototype (ICSD #55826)
        vparameters.push_back("15.489,0.47588611,0.30221447,0.204,0.152,0.31,0.337,0.322");  // 002, metallic ternary prototype (ICSD #55827)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C4_oI22_71_n_af_eh"){
        vparameters.push_back("15.1225,0.48197719,0.29965945,0.2162,0.8716,0.3059,0.3285,0.6944");  // 001, metallic ternary prototype (ICSD #156968)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C4_oI22_71_af_eh_n"){
        vparameters.push_back("14.407,0.51613799,0.2911779,0.78461,0.87104,0.6954,0.67309,0.6906");  // 001, metallic ternary prototype (ICSD #182774)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6C_oI26_71_hk_efg_a"){
        vparameters.push_back("8.613,0.97886915,0.58156275,0.655,0.665,0.336,0.757");  // 001, metallic ternary prototype (ICSD #103465)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C8_oI28_71_cg_be_2n"){
        vparameters.push_back("12.707,0.7480916,0.33721571,0.201,0.181,0.369,0.282,0.153,0.367");  // 001, metallic ternary prototype (ICSD #605107)
        vparameters.push_back("12.373,0.80029096,0.34009537,0.799,0.181,0.631,0.282,0.847,0.367");  // 002, metallic ternary prototype (ICSD #611824)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B10C_oI28_71_be_g2n_c"){
        vparameters.push_back("15.007,0.6622243,0.34248018,0.20187,0.8153,0.35777,0.71484,0.15945,0.64223");  // 001, metallic ternary prototype (ICSD #290305)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B5C3_oI28_71_gn_cn_be"){
        vparameters.push_back("12.827,0.76666407,0.34029781,0.19179,0.19303,0.3783,0.28322,0.16931,0.34743");  // 001, metallic ternary prototype (ICSD #413485)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B3C2_oI28_71_c2n_be_g"){
        vparameters.push_back("12.442,0.75952419,0.34078122,0.799,0.819,0.631,0.718,0.847,0.633");  // 001, metallic ternary prototype (ICSD #634397)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oI20_72_j_j_a"){
        vparameters.push_back("14.325,0.50820244,0.49486911,0.86413,0.70673,0.61373,0.71129");  // 001, metallic ternary prototype (ICSD #421424)
        vparameters.push_back("10.48,0.56937023,0.51545802,0.153,0.632,0.908,0.803");  // 002, ternary metal-oxide prototype (ICSD #34603)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tI10_107_a_a_ab"){
        vparameters.push_back("4.82,2.2676349,0,0.6654,0.4241,0.7494");  // 001, metallic ternary prototype (ICSD #58662) //DX20210428 - equivalent to BaNiSn3 (http://aflow.org/prototype-encyclopedia/ABC3_tI10_107_a_a_ab.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tI10_107_ab_a_a"){
        vparameters.push_back("4.2357,2.5098095,0.405,0.635,0,0.748");  // 001, metallic ternary prototype (ICSD #290390)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tI10_119_bf_a_c"){
        vparameters.push_back("4.6895,2.6930376,0.8608");  // 001, metallic ternary prototype (ICSD #249592)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_tI28_120_d_i_e"){
        vparameters.push_back("7.753,1.014962,0.856,0.648,0.625,0.417");  // 001, metallic ternary prototype (ICSD #102238)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C_tI22_121_i_ci_a"){
        vparameters.push_back("6.84,1.4307018,0.3274,0.9228,0.3237,0.2067");  // 001, metallic ternary prototype (ICSD #54354)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tP4_123_b_h_a"){
        vparameters.push_back("2.9,2.4965517,0.23");  // 001, metallic ternary prototype (ICSD #102057)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP6_123_b_i_a"){
        vparameters.push_back("4.103,1.5781136,0.3");  // 001, metallic ternary prototype (ICSD #623081)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C_tP7_123_a_ci_b"){
        vparameters.push_back("4.549,1.6322269,0.805");  // 001, metallic ternary prototype (ICSD #623944)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C_tP8_123_a_hi_b"){
        vparameters.push_back("4.35,1.8211494,0.15149,0.32935");  // 001, metallic ternary prototype (ICSD #240161)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_tP11_123_a_ehi_g"){
        vparameters.push_back("4.217,2.6013754,0.3093,0.2933,0.1189");  // 001, metallic ternary prototype (ICSD #42426)
        vparameters.push_back("4.23,2.3841608,0.3104,0.2925,0.1138");  // 002, metallic ternary prototype (ICSD #180132)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tP10_127_g_a_h"){
        vparameters.push_back("7.0659,0.48323639,0.374,0.169");  // 001, metallic ternary prototype (ICSD #54303)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP10_127_a_g_h"){
        vparameters.push_back("5.716,0.55563331,0.389,0.181");  // 001, metallic ternary prototype (ICSD #57706)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP8_129_bc_c_a"){
        vparameters.push_back("4.3365,2.2772743,0.1843,0.7494");  // 001, metallic ternary prototype (ICSD #163425)
        vparameters.push_back("4.5206,2.4320223,0.16198,0.73829");  // 002, metallic ternary prototype (ICSD #415728)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP10_129_c_ac_bc"){
        vparameters.push_back("4.401,2.2842536,0.2446,0.6257,0.8776");  // 001, metallic ternary prototype (ICSD #160053)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP10_129_c_bc_ac"){
        vparameters.push_back("4.528,2.5110424,0.75271,0.13526,0.3628");  // 001, metallic ternary prototype (ICSD #162267) //DX20210427 - equivalent to Be2CaGe2 (http://aflow.org/prototype-encyclopedia/A2BC2_tP10_129_ac_c_bc.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tP10_129_bc_c_ac"){
        vparameters.push_back("4.546,2.1240651,0.135,0.76,0.38");  // 001, metallic ternary prototype (ICSD #616562)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP14_129_ac_b_j"){
        vparameters.push_back("7.196,0.72290161,0.7755,0.43525,0.29");  // 001, metallic ternary prototype (ICSD #418547)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI8_139_b_d_a"){
        vparameters.push_back("4.137,1.7689147");  // 001, metallic ternary prototype (ICSD #54365)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI8_139_d_a_b"){
        vparameters.push_back("3.59,1.9749304");  // 001, metallic ternary prototype (ICSD #169733)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tI10_139_a_d_e"){
        vparameters.push_back("4.74,2.4978903,0.624");  // 001, metallic ternary prototype (ICSD #405)
        vparameters.push_back("3.687,2.3509628,0.8581");  // 002, ternary metal-carbo-nitride prototype (ICSD #200369) //DX20210428 - equivalent to CLi2N2 (http://aflow.org/prototype-encyclopedia/AB2C2_tI10_139_a_d_e.Li2CN2.html, part 3)
        vparameters.push_back("6.98,2.95558739255,0.208");  // 003, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI10_139_d_a_e"){
        vparameters.push_back("4.81,2.3596674,0.376");  // 001, metallic ternary prototype (ICSD #25332)
        vparameters.push_back("4.242,2.6692834,0.38820402");  // 002, metallic ternary prototype (ICSD #55789)
        vparameters.push_back("4.127,2.8279622,0.6073");  // 003, metallic ternary prototype (ICSD #57550)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI12_139_c_e_e"){
        vparameters.push_back("4.4914,3.6117024,0.66644,0.8675");  // 001, metallic ternary prototype (ICSD #182479)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI12_139_e_c_e"){
        vparameters.push_back("4.073,3.1723545,0.131,0.341");  // 001, metallic ternary prototype (ICSD #611768)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tI14_139_h_d_a"){
        vparameters.push_back("6.717,0.79082924,0.303");  // 001, metallic ternary prototype (ICSD #456)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_tI14_139_a_h_d"){
        vparameters.push_back("7.606,0.71916908,0.3");  // 001, metallic ternary prototype (ICSD #630579)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_tI16_139_bg_e_a"){
        vparameters.push_back("4.023,3.5893612,0.762,0.851");  // 001, metallic ternary prototype (ICSD #58084)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI16_139_ce_e_d"){
        vparameters.push_back("4.64,4.7327586,0.3248,0.1233");  // 001, metallic ternary prototype (ICSD #41924)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI16_139_e_ce_d"){
        vparameters.push_back("4.768,4.9496644,0.1158,0.3296");  // 001, metallic ternary prototype (ICSD #58635)
        vparameters.push_back("4.846,4.5356995,0.8717,0.6859");  // 002, metallic ternary prototype (ICSD #58638)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C3_tI20_139_e_bg_ae"){
        vparameters.push_back("4.9461,4.7316067,0.19934,0.64091,0.57874");  // 001, metallic ternary prototype (ICSD #421342)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8BC4_tI26_139_ij_a_f"){
        vparameters.push_back("8.823,0.58381503,0.6526,0.7244");  // 001, metallic ternary prototype (ICSD #57539)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B8C_tI26_139_i_fj_a"){
        vparameters.push_back("8.301,0.56161908,0.3594,0.7234");  // 001, metallic ternary prototype (ICSD #168240)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B4C_tI26_139_ij_f_a"){
        vparameters.push_back("8.86,0.5744921,0.365,0.175");  // 001, metallic ternary prototype (ICSD #57997)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_tI24_140_a_c_l"){
        vparameters.push_back("6.5562,1.721424,0.1575,0.867");  // 001, metallic ternary prototype (ICSD #172149)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI28_140_h_a_k"){
        vparameters.push_back("10.491,0.47450195,0.861,0.8027,0.0827");  // 001, metallic ternary prototype (ICSD #150145) //DX20210427 - equivalent to Sb2SiV4 (http://aflow.org/prototype-encyclopedia/A2BC4_tI28_140_h_a_k.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tI28_140_l_h_a"){
        vparameters.push_back("8.175,0.95425076,0.33913,0.36404,0.66158");  // 001, metallic ternary prototype (ICSD #182104)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_tI32_140_a_cl_h"){
        vparameters.push_back("7.7012,1.7607516,0.36124,0.65581,0.63793");  // 001, metallic ternary prototype (ICSD #107217)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_tI32_140_ah_l_c"){
        vparameters.push_back("8.565,1.4824285,0.85,0.15,0.34");  // 001, metallic ternary prototype (ICSD #616725)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_tI32_140_a_bk_h"){
        vparameters.push_back("10.586,0.48904213,0.832,0.7796,0.0724");  // 001, metallic ternary prototype (ICSD #103842)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC5_tI32_140_h_a_cl"){
        vparameters.push_back("7.9856,1.6805249,0.1336,0.6569,0.8577");  // 001, metallic ternary prototype (ICSD #156956)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_tI32_140_a_h_cl"){
        vparameters.push_back("7.4816,1.855726,0.8881,0.6658,0.8611");  // 001, metallic ternary prototype (ICSD #161658)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C3_hR14_146_ab_a2b_b"){
        vparameters.push_back("11.53099,0.4552076,0.944,0.46,0.352,0.533,0.121,0.408,0.149,0.972,0.617,0.864,0.026,0.629,0.461,0.903");  // 001, metallic ternary prototype (ICSD #103911)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_156_3a_b2c_2bc"){
        vparameters.push_back("4.94,2.2064777,0.334,0.668,0,0.454,0.146,0.853,0.214,0.788,0.5");  // 001, metallic ternary prototype (ICSD #58911)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC4_hR10_160_5a_a_4a"){
        vparameters.push_back("4.7099904,6.7346082,0.94,0.862,0.409,0.061,0.145,0.501,0.7787,0.6863,0.2327,0.3215");  // 001, metallic ternary prototype (ICSD #12142)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hR4_166_c_b_a"){
        vparameters.push_back("4.2709972,3.7286311,0.1888");  // 001, metallic ternary prototype (ICSD #102938)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_hR7_166_a_c_2c"){
        vparameters.push_back("4.2259957,6.2863209,0.90091,0.79965,0.61204");  // 001, metallic ternary prototype (ICSD #150127)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C6_hR11_166_a_2c_h"){
        vparameters.push_back("5.9937387,3.892206,0.66786054,0.87208523,0.39951361,0.88682193");  // 001, metallic ternary prototype (ICSD #166193)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hR12_166_h_bc_ac"){
        vparameters.push_back("5.4710243,4.6349773,0.335,0.142,0.593,0.054");  // 001, metallic ternary prototype (ICSD #604213)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B3C2_hR12_166_bh_ac_c"){
        vparameters.push_back("5.6140311,4.6050915,0.1426,0.335,0.5787,0.0745");  // 001, metallic ternary prototype (ICSD #57538)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C9_hR12_166_a_c_bch"){
        vparameters.push_back("4.9234811,4.8473628,0.146,0.334,0.5872,0.0827");  // 001, metallic ternary prototype (ICSD #55614)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2C_hR12_166_eh_c_a"){
        vparameters.push_back("5.1579625,3.6593664,0.222,0.278,0.777");  // 001, metallic ternary prototype (ICSD #646978)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR12_166_bc_h_ac"){
        vparameters.push_back("5.6320019,4.8020254,0.333,0.14,0.586,0.074");  // 001, metallic ternary prototype (ICSD #604688)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C3_hR12_166_c_bh_ac"){
        vparameters.push_back("5.6139959,4.7716445,0.6644,0.8571,0.419,0.908");  // 001, metallic ternary prototype (ICSD #104173)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B4C2_hR13_166_ah_2c_c"){
        vparameters.push_back("6.8009851,4.2773178,0.3152,0.1878,0.5493,0.2133,0.8046");  // 001, metallic ternary prototype (ICSD #58581)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C2_hR14_166_3c_2c_abc"){
        vparameters.push_back("4.7683821,9.8351827,0.93576,0.79728,0.64381,0.87437,0.70884,0.57697");  // 001, metallic ternary prototype (ICSD #710044)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_b_a_b"){
        vparameters.push_back("4.7199,1.6142715,0.75,0.0593,0.473");  // 001, metallic ternary prototype (ICSD #54997)
        vparameters.push_back("4.576,1.6121066,0,0.382,0");  // 002, metallic ternary prototype (ICSD #100115)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_186_a_b_b"){
        vparameters.push_back("4.5858,1.7114789,0.25,0.523,0.978");  // 001, metallic ternary prototype (ICSD #156394)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_hP12_186_ac_b_b"){
        vparameters.push_back("4.9708,1.5973284,0,0.9386,0.5644,0.83026667,0.7503");  // 001, metallic ternary prototype (ICSD #424277)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_hP16_186_a2b_b_2a2b"){
        vparameters.push_back("4.4831,4.5854431,0.16473,0.02937,0.29581,0.98426,0.84576,0.711,0.12513,0.40817");  // 001, metallic ternary prototype (ICSD #412207)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_187_i_h_ab"){
        vparameters.push_back("4.3442,1.5409512,0.7853,0.7161");  // 001, metallic ternary prototype (ICSD #156264)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_hP6_187_hi_b_a"){
        vparameters.push_back("4.3168,1.5918273,0.30839,0.22249");  // 001, metallic ternary prototype (ICSD #98666)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_187_dh_ai_eg"){
        vparameters.push_back("4.933,2.2278532,0.6533,0.8342,0.708");  // 001, metallic ternary prototype (ICSD #409533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B5C_hP11_187_hk_cgi_a"){
        vparameters.push_back("6.0692,1.8322019,0.32828,0.322,0.2,0.5122");  // 001, metallic ternary prototype (ICSD #416337)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C_hP11_187_ck_jk_a"){
        vparameters.push_back("8.0581,0.55113488,0.53424,0.79441,0.19539");  // 001, metallic ternary prototype (ICSD #249520)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C3_hP8_189_f_c_g"){
        vparameters.push_back("7.2873,0.42088291,0.76243633,0.41574233");  // 001, metallic ternary prototype (ICSD #410967)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_g_f_ad"){
        vparameters.push_back("7.579,0.51801029,0.60163,0.26746");  // 001, metallic ternary prototype (ICSD #51845) //DX20210428 - equivalent to AlNiZr (http://aflow.org/prototype-encyclopedia/ABC_hP9_189_g_ad_f.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C6_hP9_189_b_c_fg"){
        vparameters.push_back("7.801,0.43276503,0.249,0.592");  // 001, metallic ternary prototype (ICSD #20876)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3C_hP9_189_cf_g_b"){
        vparameters.push_back("6.759,0.61236869,0.24,0.578");  // 001, metallic ternary prototype (ICSD #634575)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_g_bc_f"){
        vparameters.push_back("7.717,0.49941687,0.597,0.26");  // 001, metallic ternary prototype (ICSD #54315)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_hP9_189_c_b_fg"){
        vparameters.push_back("8.213,0.50651406,0.235,0.601");  // 001, metallic ternary prototype (ICSD #96253)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_bc_g"){
        vparameters.push_back("7.637,0.61332984,0.75,0.4277");  // 001, metallic ternary prototype (ICSD #54344)
        vparameters.push_back("6.798,0.46587232,0.273,0.615");  // 002, metallic ternary prototype (ICSD #103885)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_f_g_bc"){
        vparameters.push_back("7.0015,0.54982504,0.793,0.446");  // 001, metallic ternary prototype (ICSD #107416)
        vparameters.push_back("7.0941,0.50138848,0.76,0.6667");  // 002, metallic ternary prototype (ICSD #608141)
        vparameters.push_back("7.016,0.5791049,0.24,0.3333");  // 003, metallic ternary prototype (ICSD #608770)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_bc_f_g"){
        vparameters.push_back("7.45,0.53288591,0.75506667,0.51506667");  // 001, metallic ternary prototype (ICSD #628137)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C_hP6_191_g_c_a"){
        vparameters.push_back("5.565,0.81239892");  // 001, metallic ternary prototype (ICSD #57329)
        vparameters.push_back("8.813,0.46102349");  // 002, metallic ternary prototype (ICSD #102452)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP6_191_c_a_g"){
        vparameters.push_back("7.44,0.52956989");  // 001, metallic ternary prototype (ICSD #108369)
        vparameters.push_back("5.581,0.70757929");  // 002, metallic ternary prototype (ICSD #658142)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9BC2_hP12_191_fm_a_c"){
        vparameters.push_back("8.04,0.48383085,0.212");  // 001, metallic ternary prototype (ICSD #57518)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC6_hP13_191_i_a_cde"){
        vparameters.push_back("5.48,1.6394161,0.6683,0.7529");  // 001, metallic ternary prototype (ICSD #54273)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6C_hP13_191_i_cde_a"){
        vparameters.push_back("5.5301,1.6316522,0.168,0.746");  // 001, metallic ternary prototype (ICSD #54277)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C2_hP16_193_g_g_d"){
        vparameters.push_back("7.9396,0.68974508,0.60359025,0.24573816");  // 001, metallic ternary prototype (ICSD #103733)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_d_a_c"){
        vparameters.push_back("4.713,1.8949714");  // 001, metallic ternary prototype (ICSD #54319)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_c_a_d"){
        vparameters.push_back("4.553,1.7665276");  // 001, metallic ternary prototype (ICSD #57018)
        vparameters.push_back("4.5813,1.6472617");  // 002, metallic ternary prototype (ICSD #600994)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_c_d"){
        vparameters.push_back("4.801,2.0249948");  // 001, metallic ternary prototype (ICSD #106315)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP6_194_a_d_c"){
        vparameters.push_back("4.909,1.500713");  // 001, metallic ternary prototype (ICSD #106345)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_b_f_c"){
        vparameters.push_back("4.267,2.0180455,0.58");  // 001, metallic ternary prototype (ICSD #58139)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_a_f_c"){
        vparameters.push_back("4.539,1.9649703,0.592");  // 001, metallic ternary prototype (ICSD #59505)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hP8_194_f_a_c"){
        vparameters.push_back("4.303,1.7748083,0.583");  // 001, metallic ternary prototype (ICSD #150602)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_hP12_194_h_a_f"){
        vparameters.push_back("5.27,1.626186,0.063,0.833");  // 001, metallic ternary prototype (ICSD #58159)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_194_f_f_ab"){
        vparameters.push_back("4.552,4.1001757,0.1617,0.6162");  // 001, metallic ternary prototype (ICSD #66003)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP12_194_ab_f_f"){
        vparameters.push_back("4.4951,3.5496874,0.11481,0.65681");  // 001, metallic ternary prototype (ICSD #152621)
        vparameters.push_back("4.588,3.7149085,0.6154,0.1507");  // 002, metallic ternary prototype (ICSD #152623)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C2_hP16_194_h_af_f"){
        vparameters.push_back("5.9892,2.4310425,0.12318,0.612,0.84546667");  // 001, metallic ternary prototype (ICSD #415931)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cP16_215_e_bce_ad"){
        vparameters.push_back("6.315,0.7511,0.2331");  // 001, metallic ternary prototype (ICSD #180119)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_cF16_216_a_d_bc"){
        vparameters.push_back("6.35");  // metallic ternary prototype (ICSD #57330) //DX20210428 - equivalent to CuHg2Ti inverse Heulser (http://aflow.org/prototype-encyclopedia/AB2C_cF16_216_b_ad_c.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_cF24_216_a_e_c"){
        vparameters.push_back("6.86,0.625");  // 001, metallic ternary prototype (ICSD #59462)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C4_cI22_217_b_c_c"){
        vparameters.push_back("7.41,0.2031,0.366");  // 001, metallic ternary prototype (ICSD #58899)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_cP15_221_ag_c_d"){
        vparameters.push_back("6.408,0.29015");  // 001, metallic ternary prototype (ICSD #412144)
      }

      // -------------------------------------------------------------------------
      // metal-boride prototypes (from DX) //DX20210104
      // ---------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_mC10_12_j_a"){
        vparameters.push_back("5.476,0.9800767,0.84556245,147.4208,0.5019,0.8429,0.8071");  // 001, binary metal-boride prototype (ICSD #15079)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_mC28_15_ef_2f"){
        vparameters.push_back("6.4282,0.75907719,1.38538,121.3067,0.8142,0.9572,0.827,0.4396,0.6893,0.0019,0.4835,0.334,0.6805,0.2864");  // 001, binary metal-boride prototype (ICSD #24308)
        vparameters.push_back("6.43,0.7592535,1.3852411,121.3305,0.67,0.96,0.66,0.44,0.69,0,0.483,0.335,0.682,0.287");  // 002, binary metal-boride prototype (ICSD #150561)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oP10_58_2g_a"){
        vparameters.push_back("4.7452,1.1540926,0.60402091,0.8357,0.3669,0.7764,0.6791");  // 001, binary metal-boride prototype (ICSD #186851)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP6_59_e_a"){
        vparameters.push_back("2.8668,1.6200293,1.410946,0.8477,0.4373,0.3691");  // 001, binary metal-boride prototype (ICSD #31871) //DX20210428 - equivalent to B2Ru1 (http://aflow.org/prototype-encyclopedia/A2B_oP6_59_f_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_62_d_c"){
        vparameters.push_back("4.8159,0.99811043,0.7765319,0.51985,0.62385,0.66113,0.4307,0.10703");  // 001, binary metal-boride prototype (ICSD #425310)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC8_63_a_c"){
        vparameters.push_back("3.3,1.7212121,1.3151515,0.333");  // 001, binary metal-boride prototype (ICSD #150732)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_63_c_cf"){
        vparameters.push_back("2.89,3.2224913,2.5114187,0.744,0.4262,0.8655,0.938");  // 001, binary metal-boride prototype (ICSD #43662) //DX20210428 - equivalent to Re3B (http://aflow.org/prototype-encyclopedia/AB3_oC16_63_c_cf.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oC20_63_3c_2c"){
        vparameters.push_back("3.0599,6.0227458,0.97516259,0.97649,0.88221,0.16862,0.57065,0.29502");  // 001, binary metal-boride prototype (ICSD #79258)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B5_oC22_65_i2j_aij"){
        vparameters.push_back("3.1385,7.2015294,1.0481122,0.7266,0.61446,0.9168,0.5397,0.80566");  // 001, binary metal-boride prototype (ICSD #68538)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oI10_71_n_a"){
        vparameters.push_back("5.4773,0.86613843,0.52323225,0.3455,0.8249");  // 001, binary metal-boride prototype (ICSD #24353)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oI14_71_ef_af"){
        vparameters.push_back("14.82,0.22300945,0.21167341,0.556,0.124,0.32");  // 001, binary metal-boride prototype (ICSD #76631)
        vparameters.push_back("13.4916,0.25517359,0.22686709,0.43285384,0.87121302,0.6865067");  // 002, binary metal-boride prototype (ICSD #614732)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B_tI24_119_aegi_d"){
        vparameters.push_back("4.404,1.7502271,0.3298,0.2922,0.2922,0.8339");  // 001, binary metal-boride prototype (ICSD #164844)
        vparameters.push_back("4.984,1.168138,0.7154,0.6847,0.3188,0.867");  // 002, binary metal-boride prototype (ICSD #164845)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI12_121_ab_i"){
        vparameters.push_back("5.098,0.83169871,0.6667,0.7");  // 001, binary metal-boride prototype (ICSD #16809) //DX20210428 - equivalent to Fe2B (C17, obsolete, http://aflow.org/prototype-encyclopedia/AB2_tI12_121_ab_i.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_tI26_139_jm_a"){
        vparameters.push_back("5.2347,1.4056775,0.837,0.832,0.668");  // 001, binary metal-boride prototype (ICSD #615424)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI12_140_b_h"){
        vparameters.push_back("5.11,0.83150685,0.245");  // 001, binary metal-boride prototype (ICSD #160789)
        vparameters.push_back("5.11,0.83150685,0.315");  // 002, binary metal-boride prototype (ICSD #160790)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI12_140_c_h"){
        vparameters.push_back("5.647,0.83956083,0.3333");  // 001, binary metal-boride prototype (ICSD #189385)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5_hR9_160_ab_2ab"){
        vparameters.push_back("6.9720729,1.2247449,0.075,0.825,0.325,0.175,0.825,0.675,0.325");  // 001, binary metal-boride prototype (ICSD #23041)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B_hR12_160_2a3b_a"){
        vparameters.push_back("4.9820288,2.2326381,0.0127,0.7936,0.4656,0.3767,0.0356,0.0203,0.379,0.7862,0.4275");  // 001, binary metal-boride prototype (ICSD #164842)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hR6_166_2c_c"){
        vparameters.push_back("3.0136115,6.9481798,0.3323,0.18184,0.07569");  // 001, binary metal-boride prototype (ICSD #39554)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR8_166_f_c"){
        vparameters.push_back("5.2239951,1.7923052,0.1649,0.6651");  // 001, binary metal-boride prototype (ICSD #167734)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP3_191_d_a"){
        vparameters.push_back("3.1593,1.4315513");  // 001, binary metal-boride prototype (ICSD #169458)
        vparameters.push_back("3.163,1.2611445");  // 002, binary metal-boride prototype (ICSD #186764)
        vparameters.push_back("3.1,0.74096774");  // 003, binary metal-boride prototype (ICSD #614887)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_194_f_c"){
        vparameters.push_back("2.9,2.5786207,0.548");  // 001, binary metal-boride prototype (ICSD #23871)
        vparameters.push_back("2.583,3.8946961,0.549");  // 002, binary metal-carbide prototype (ICSD #184660)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_194_bc_f"){
        vparameters.push_back("3.057,3.6113837,0.626");  // 001, binary metal-boride prototype (ICSD #2)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_194_af_c"){
        vparameters.push_back("2.9,2.5775862,0.55");  // 001, binary metal-boride prototype (ICSD #24361) //DX20210428 - equivalent to B3Re (http://aflow.org/prototype-encyclopedia/A3B_hP8_194_af_c.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_hP10_194_2f_c"){
        vparameters.push_back("2.951,3.7217892,0.0449,0.6114");  // 001, binary metal-boride prototype (ICSD #167735)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP10_194_bf_f"){
        vparameters.push_back("3.0797,8.2163847,0.31026779,0.86352537");  // 001, binary metal-boride prototype (ICSD #108082)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP12_194_bcf_f"){
        vparameters.push_back("2.9831,4.6525427,0.4757,0.63759");  // 001, binary metal-boride prototype (ICSD #23716)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_hP14_194_fg_e"){
        vparameters.push_back("3.015,4.907131,0.3981,0.3023");  // 001, binary metal-boride prototype (ICSD #167733)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B_cF48_216_aeg_b"){
        vparameters.push_back("6.81,0.8481,0.442");  // 001, binary metal-boride prototype (ICSD #164841)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B_cF48_216_aef_d"){
        vparameters.push_back("6.703,0.1534,0.6772");  // 001, binary metal-boride prototype (ICSD #164843)
      }
      // ---------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_aP7_2_2i_i_a"){
        vparameters.push_back("3.1862,2.0826376,0.93500094,77.1999,90.3494,79.8085,0.63794635,0.7241073,0.13794635,0.56918124,0.86163752,0.56918124,0.80787091,0.38425819,0.30787091");  // 001, ternary metal-boride prototype (ICSD #44175)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C3_mC20_5_b2c_c_ac"){
        vparameters.push_back("10.4164,0.38525786,0.8516282,131.9992,0,0.205,0.911,0.492,0.88,0.163,0.205,0.308,0.3987,0.213,0.2398,0.2683,0.7352,0.4025");  // 001, ternary metal-boride prototype (ICSD #59227)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_mP10_11_3e_e_e"){
        vparameters.push_back("6.79,0.46410898,0.79469809,101.6,0.525,0.584,0.008,0.59,0.004,0.915,0.7091,0.1764,0.2217,0.3097");  // 001, ternary metal-boride prototype (ICSD #65932)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_mC12_12_g_a_df"){
        vparameters.push_back("5.409,1.7339619,0.57330375,91.2,0.333");  // 001, ternary metal-boride prototype (ICSD #44230)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C5_mC22_12_2i_i_aj"){
        vparameters.push_back("9.786,0.5270795,0.82546495,136.63,0.119,0.341,0.041,0.639,0.4228,0.1487,0.7964,0.751,0.3707");  // 001, ternary metal-boride prototype (ICSD #63501)
        vparameters.push_back("9.7745,0.55030948,0.84464679,137.196,0.119,0.341,0.041,0.639,0.4167,0.1649,0.784,0.75,0.362");  // 002, ternary metal-boride prototype (ICSD #170618)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_mC20_15_f_f_e"){
        vparameters.push_back("8.389,0.62522351,0.82963404,127.6585,0.662,0.101,0.924,0.022,0.321,0.635,0.199");  // 001, ternary metal-boride prototype (ICSD #57003)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oF48_43_b_b_b"){
        vparameters.push_back("9.963,0.87604135,0.65070762,0.485,0.307,0.478,0.3641,0.8374,0.6596,0.13652,0.4853,0");  // 001, ternary metal-boride prototype (ICSD #75029)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP6_47_d_rt_a"){
        vparameters.push_back("3.1346,1.1159319,1.7320551,0.3333,0.1667");  // 001, ternary metal-boride prototype (ICSD #181368)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C3_oP8_47_a_qr_dt"){
        vparameters.push_back("2.952,1.0125339,2.7408537,0.288,0.394,0.188");  // 001, ternary metal-boride prototype (ICSD #20082)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_oP16_47_qrs_qrst_ad"){
        vparameters.push_back("2.7998,1.0070719,8.2780913,0.1224,0.30476,0.235,0.41598,0.6224,0.80476,0.91598");  // 001, ternary metal-boride prototype (ICSD #71640)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3_oP14_55_h_h_ag"){
        vparameters.push_back("5.773,1.6284428,0.48865408,0.6656,0.7789,0.938,0.832,0.843,0.574");  // 001, ternary metal-boride prototype (ICSD #417442)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP16_62_2c_c_c"){
        vparameters.push_back("6.057,0.51626218,1.356282,0.758,0.596,0.83,0.192,0.005,0.38,0.865,0.87");  // 001, ternary metal-boride prototype (ICSD #41893)
        vparameters.push_back("5.636,0.59811923,1.9201561,0.844,0.248,0.247,0.994,0.918,0.852,0.582,0.102");  // 002, ternary metal-nitride prototype (ICSD #617606)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP16_62_d_c_c"){
        vparameters.push_back("5.8948,0.90983579,1.0952195,0.989,0.335,0.818,0.818,0.642,0.084,0.536");  // 001, ternary metal-boride prototype (ICSD #612853)
        vparameters.push_back("8.3113,0.62342834,1.1352255,0.5088,0.816,0.1472,0.0799,0.77,0.134,0.11");  // 002, ternary metal-nitride prototype (ICSD #189822)
        vparameters.push_back("8.578,0.59552343,1.1488109,0.5069,0.8213,0.1584,0.0802,0.787,0.133,0.093");  // 003, ternary metal-nitride prototype (ICSD #189823)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_oC20_63_3c_c_c"){
        vparameters.push_back("3.04,5.7730263,0.98026316,0.674,0.475,0.347,0.795,0.072");  // 001, ternary metal-boride prototype (ICSD #44188)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_oC32_63_2cf_a_cf"){
        vparameters.push_back("2.9128,4.0713403,3.3500412,0.845,0.123,0.4851,0.656,0.891,0.7002,0.3927");  // 001, ternary metal-boride prototype (ICSD #91243)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_oC10_65_a_i_j"){
        vparameters.push_back("2.937,3.7691522,1.0115764,0.22,0.352");  // 001, ternary metal-boride prototype (ICSD #20083)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_oC12_65_i_a_cf"){
        vparameters.push_back("5.5452,1.7372142,0.55383034,0.33330607");  // 001, ternary metal-boride prototype (ICSD #40777)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3_oC14_65_j_j_ai"){
        vparameters.push_back("2.967,5.7418268,0.99932592,0.17855,0.4174,0.2768");  // 001, ternary metal-boride prototype (ICSD #43843)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C3_oC22_65_gp_j_ah"){
        vparameters.push_back("7.701,1.1208934,0.44955201,0.674,0.7754,0.2728,0.766,0.157");  // 001, ternary metal-boride prototype (ICSD #65737)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oC24_65_i_bc_jkm"){
        vparameters.push_back("5.0411,1.7389459,1.374819,0.16661035,0.16692044,0.79646249,0.79638046");  // 001, ternary metal-boride prototype (ICSD #20882)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_65_h_gip_j"){
        vparameters.push_back("7.308,1.2746305,0.47865353,0.376,0.8184,0.908,0.69941,0.7768,0.8391");  // 001, ternary metal-boride prototype (ICSD #167533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C5_oF44_69_o_g_acf"){
        vparameters.push_back("11.22,0.86693405,0.48663102,0.3576,0.3701,0.3396");  // 001, ternary metal-boride prototype (ICSD #66768)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oF40_70_g_e_a"){
        vparameters.push_back("10.672,0.56118816,0.88071589,0.7551,0.453");  // 001, ternary metal-boride prototype (ICSD #8153)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oF40_70_f_a_e"){
        vparameters.push_back("10.0613,0.90966376,0.6569827,0.76165914,0.80936151");  // 001, ternary metal-boride prototype (ICSD #601369)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oF48_70_e_f_c"){
        vparameters.push_back("9.154,1.1433253,0.66812322,0.802,0.49266");  // 001, ternary metal-boride prototype (ICSD #75027)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oI10_71_h_a_f"){
        vparameters.push_back("7.1344,0.64064252,0.44652669,0.20187042,0.70789234");  // 001, ternary metal-boride prototype (ICSD #16776)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_oI14_71_ef_f_a"){
        vparameters.push_back("12.7,0.24629921,0.23496063,0.556,0.125,0.32");  // 001, ternary metal-boride prototype (ICSD #44293)
        vparameters.push_back("13.528,0.23945151,0.22450473,0.43188933,0.86540104,0.68747521");  // 002, ternary metal-boride prototype (ICSD #614782)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC2_oI18_71_hn_a_f"){
        vparameters.push_back("8.318,0.78865112,0.37292618,0.848,0.635,0.67,0.29");  // 001, ternary metal-boride prototype (ICSD #16203)
        vparameters.push_back("8.341,0.78671622,0.37225752,0.84567,0.6438,0.3035,0.2373");  // 002, ternary metal-boride prototype (ICSD #81544)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tP5_99_a_b_ac"){
        vparameters.push_back("4.0031,1.2674178,0.688,0.1174,0,0.5132");  // 001, ternary metal-boride prototype (ICSD #95049)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI10_139_e_a_d"){
        vparameters.push_back("3.915,2.8735632,0.319");  // 001, ternary metal-boride prototype (ICSD #8155)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_tI10_139_d_e_a"){
        vparameters.push_back("2.9631,4.259323,0.62653");  // 001, ternary metal-boride prototype (ICSD #71642)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_tI10_139_e_d_a"){
        vparameters.push_back("3.5369,2.6189318,0.65877134");  // 001, ternary metal-boride prototype (ICSD #612907) //DX20210428 - equivalent to Cr2Si2Th (http://aflow.org/prototype-encyclopedia/A2B2C_tI10_139_d_e_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C7_hR12_166_c_bc_ah"){
        vparameters.push_back("5.1569859,4.1560984,0.831,0.6449,0.06106,0.5723");  // 001, ternary metal-boride prototype (ICSD #36505)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B3C4_hR13_166_g_ac_2c"){
        vparameters.push_back("5.481,4.4389657,0.889,0.748,0.584,0.835");  // 001, ternary metal-boride prototype (ICSD #614152)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_hP16_176_bh_c_h"){
        vparameters.push_back("7.5482,0.46199094,0.95606667,0.39513333,0.74496667,0.073453333");  // 001, ternary metal-boride prototype (ICSD #1518)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP12_180_d_c_i"){
        vparameters.push_back("5.4898,1.4364822,0.15128");  // 001, ternary metal-boride prototype (ICSD #90403)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_hP9_189_ac_f_g"){
        vparameters.push_back("6.015,0.53566085,0.74406667,0.40006667");  // 001, ternary metal-boride prototype (ICSD #20298)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C5_hP13_189_ab_i_df"){
        vparameters.push_back("6.548,0.85113012,0.5982,0.2551,0.2542");  // 001, ternary metal-boride prototype (ICSD #77352)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP15_189_bc_f_fk"){
        vparameters.push_back("9.236,0.29969684,0.66866667,0.17796667,0.83806667,0.32113333");  // 001, ternary metal-boride prototype (ICSD #68091)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_191_c_ab_i"){
        vparameters.push_back("5.547,0.99224806,0.7422");  // 001, ternary metal-boride prototype (ICSD #68092)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hP12_191_d_ci_ab"){
        vparameters.push_back("4.4977,1.5434555,0.29");  // 001, ternary metal-boride prototype (ICSD #16515)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_hP12_191_c_ab_di"){
        vparameters.push_back("5.066,1.3509672,0.2054");  // 001, ternary metal-boride prototype (ICSD #36504)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C_cF40_225_c_bd_a"){
        vparameters.push_back("8.23");  // ternary metal-boride prototype (ICSD #260147)
      }

      // -------------------------------------------------------------------------
      // metal-oxide prototypes (from DX) //DX20210106
      // ---------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mP12_4_2a_4a"){
        vparameters.push_back("5.6004,0.86475966,0.98671523,119.4472,0.965,0.5,0.733,0.465,0,0.233,0.84,0.82,0.45,0.12,0.22,0.01,0.34,0.72,0.95,0.62,0.3,0.51");  // 001, binary metal-oxide prototype (ICSD #36263)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_mC10_5_2c_a"){
        vparameters.push_back("8.66,0.52193995,0.88742494,146.891,0,0.8,0.79,0.93,0.2,0.21,0.31");  // 001, binary metal-oxide prototype (ICSD #24672)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5_mC14_5_c_a2c"){
        vparameters.push_back("5.2193,0.9004081,1.2518729,120.6649,0.282,0.0553,0,0.7654,0.278,0.326,0.437,0.73,0.827,0.806");  // 001, binary metal-oxide prototype (ICSD #51176)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC6_8_a_2a"){
        vparameters.push_back("4.8339,0.58416599,2.3231138,107.0815,0,0,0.721,0.084,0.26,0.896");  // 001, binary metal-oxide prototype (ICSD #95440)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC8_9_a_a"){
        vparameters.push_back("4.6927,0.73056023,1.3546999,127.1683,0.75,0.5033,0,0.7706,0.8323,0.2602");  // 001, binary metal-oxide prototype (ICSD #69757)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP6_10_mn_ce"){
        vparameters.push_back("4.617,0.6278969,0.97595841,91.79,0.294,0.239,0.796,0.274");  // 001, binary metal-oxide prototype (ICSD #89471)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_10_2mo_im"){
        vparameters.push_back("9.06,0.6401766,0.49908389,91.85,0.7189,0.3975,0.2284,0.098,0.7862,0.7688,0.4689,0.854,0.2474,0.7135");  // 001, binary metal-oxide prototype (ICSD #1503)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mP8_11_e_3e"){
        vparameters.push_back("7.095,0.51966173,0.55729387,103.75,0.2966,0.3286,0.6309,0.5599,0.324,0.9076,0.0547,0.2709");  // 001, binary metal-oxide prototype (ICSD #80577)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_mP10_11_3e_2e"){
        vparameters.push_back("7.064,0.79657418,0.54827293,99.99,0.2278,0.2131,0.6194,0.2486,0.9306,0.65,0.9167,0.1444,0.3833,0.6722");  // 001, binary metal-oxide prototype (ICSD #36243)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mP12_11_4e_2e"){
        vparameters.push_back("12.0123,0.31144743,0.5431849,104.236,0.868,0.872,0.736,0.389,0.94,0.311,0.638,0.931,0.803,0.095,0.901,0.606");  // 001, binary metal-oxide prototype (ICSD #657748)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_mP14_11_3e_4e"){
        vparameters.push_back("10.88,0.31709559,0.69944853,109.58,0.6158,0.6443,0.1769,0.4707,0.5967,0.2381,0.8673,0.7738,0.8538,0.3059,0.3665,0.0699,0.3914,0.5744");  // 001, binary metal-oxide prototype (ICSD #23478)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_mP14_11_5e_2e"){
        vparameters.push_back("7.114,0.5020804,0.88341299,90.069,0.8241,0.4374,0.0307,0.8536,0.5086,0.6752,0.6988,0.0424,0.1992,0.2627,0.9002,0.1897,0.7227,0.7443");  // 001, binary metal-oxide prototype (ICSD #59960)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC12_12_i_2i"){
        vparameters.push_back("6.6495,0.60415069,0.98013384,81.89,0.20088,0.76934,0.5771,0.5709,0.1544,0.1381");  // 001, binary metal-oxide prototype (ICSD #62228)
        vparameters.push_back("6.082,0.57875699,2.1769155,108.2,0.1872,0.7549,0.5798,0.8696,0.7548,0.6305");  // 002, binary metal-oxide prototype (ICSD #77699)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC18_12_ai_3i"){
        vparameters.push_back("13.7,0.20927007,0.32554745,89.5,0.339,0.506,0.5535,0.768,0.1095,0.26,0.217,0.726");  // 001, binary metal-oxide prototype (ICSD #150462)
        vparameters.push_back("14.73541,0.247159732915,0.636222541483,110.41211,0.36025,0.7087,0.3503,-0.0096,0.703,0.6208,0.0039,0.3098");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_mC20_12_di_a3i"){
        vparameters.push_back("9.507,0.30871989,0.8094036,89.16,0.788,0.637,0.867,0.388,0.583,0.251,0.295,0.115");  // 001, binary metal-oxide prototype (ICSD #77706)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7_mC20_12_ai_d3i"){
        vparameters.push_back("9.5536,0.30547647,0.81248953,89.68,0.7143,0.1429,0.3572,0.0714,0.2143,0.6428,0.0714,0.2143");  // 001, binary metal-oxide prototype (ICSD #163494)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8_mC22_12_ai_2ij"){
        vparameters.push_back("8.4067,0.62635755,0.71921206,88.97,0.5528,0.2402,0.2007,0.9116,0.2971,0.523,0.5495,0.2567,0.1917");  // 001, binary metal-oxide prototype (ICSD #155847)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC24_12_4i_2i"){
        vparameters.push_back("12.03,0.30698254,0.99438903,148.9606,0.128,0.991,0.635,0.373,0.661,0.595,0.087,0.729,0.92,0.725,0.398,0.3");  // 001, binary metal-oxide prototype (ICSD #199)
        vparameters.push_back("12.16,0.30756579,0.9840625,148.7048,0.38,0.01,0.1,0.34,0.19,0.63,0.57,0.72,0.6,0.28,0.86,0.26");  // 002, binary metal-oxide prototype (ICSD #57154)
        vparameters.push_back("12.093,0.30613578,0.98613247,148.9388,0.6399,0,0.1098,0.3436,0.2064,0.6496,0.5716,0.6928,0.4204,0.7214,0.915,0.3145");  // 003, binary metal-oxide prototype (ICSD #73855)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC24_12_2ij_gi"){
        vparameters.push_back("9.0664,0.63939381,0.49915071,88.12,0.2811,0.6031,0.2089,0.9,0.7987,0.7686,0.5312,0.3518,0.2525,0.2942");  // 001, binary metal-oxide prototype (ICSD #34416)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B8_mC26_12_ahi_2ij"){
        vparameters.push_back("10.347,0.55320383,0.95300087,152.3489,0.238,0.568,0.838,0.688,0.595,0.186,0.575,0.704,0.266,0.106");  // 001, binary metal-oxide prototype (ICSD #16956)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_mP14_13_2g_ag"){
        vparameters.push_back("8,0.6875,0.6,95.9,0.77,0.68,0.5,0.1,0.7,0.23,0.67,0.5,0.29");  // 001, binary metal-oxide prototype (ICSD #174299)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP8_14_ad_e"){
        vparameters.push_back("5.852,0.59432673,0.93899522,107.5,0.705,0.35,0.77");  // 001, binary metal-oxide prototype (ICSD #43741)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_mP14_14_ae_2e"){
        vparameters.push_back("5.809,1.5851093,0.97729385,143.7156,0.7151,0.642,0.7632,0.5121,0.5499,0.2857,0.9634,0.8119,0.7878");  // 001, binary metal-oxide prototype (ICSD #59225)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC12_15_c_f"){
        vparameters.push_back("7.88,0.51218274,0.96200508,62.0091,0.44,0.3,0.78");  // 001, binary metal-oxide prototype (ICSD #38230)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mC24_15_ae_2f"){
        vparameters.push_back("12.3668,0.41384999,0.96273086,153.57,0.482,0.2233,0.687,0.652,0.6297,0.334,0.973");  // 001, binary metal-oxide prototype (ICSD #79500)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_mC32_15_e2f_cf"){
        vparameters.push_back("9.98,0.50400802,0.98597194,138.8,0.9,0.206,0.1,0.148,0.111,0.6,0.939,0.162,0.75,0.786");  // 001, binary metal-oxide prototype (ICSD #15899)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_mC32_15_e2f_af"){
        vparameters.push_back("9.846,0.51054235,0.71186269,70.464,0.80642,0.80743,0.1525,0.03331,0.919,0.65599,0.13796,0.62965,0.00159,0.58985");  // 001, binary metal-oxide prototype (ICSD #32587)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B3_oC22_21_achl_bg"){
        vparameters.push_back("6.386,1.6332603,0.59505168,0.824,0.802,0.343,0.152,0.85");  // 001, binary metal-oxide prototype (ICSD #73719)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_29_a_a"){
        vparameters.push_back("5.48,0.8649635,1.0729927,0.85,0.3583,0.1111,0.98611,0.25,0.7694");  // 001, binary metal-oxide prototype (ICSD #36250)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP4_31_a_a"){
        vparameters.push_back("3.82,0.94502618,1.1256545,0.25,0.25,0.5,0.8");  // 001, binary metal-oxide prototype (ICSD #20624)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_oP14_31_a2b_b"){
        vparameters.push_back("11.503,0.37981396,0.30922368,0.92,0.89,0.852,0.55,0.92,0.8,0.97,0.46,0.852,0.903,0");  // 001, binary metal-oxide prototype (ICSD #15984) //DX20210428 - equivalent to O5V2 (shcherbinaite (obsolete), D8_{7}, http://aflow.org/prototype-encyclopedia/A5B2_oP14_31_a2b_b.html, part 3)
        vparameters.push_back("11.48,0.37979094,0.30966899,0.83,0.11,0.845,0.55,0.92,0.805,0.96,0.46,0.854,0.905,0");  // 002, binary metal-oxide prototype (ICSD #29140)
        vparameters.push_back("11.555,0.3782778,0.30783211,0.003,0.002,0.8532,0.534,0.032,0.8187,0.0056,0.5069,0.8513,0.8914,0");  // 003, binary metal-oxide prototype (ICSD #41030)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oC24_35_abef_de"){
        vparameters.push_back("6.6941,1.6509912,1.9077845,0.14824608,0.056987889,0.77696186,0.062716618,0.28327908,0.14033669,0.27745706,0.01111837,0.79079507,0.86091387,0.98198067");  // 001, binary metal-oxide prototype (ICSD #97008)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oC16_36_b_2a"){
        vparameters.push_back("5,1.144,2.224,0,0,0.559,0.254,0.32,0.25,0.38");  // 001, binary metal-oxide prototype (ICSD #60619)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B_oC18_38_2cde_a"){
        vparameters.push_back("6.02,0.90863787,1.551495,0,0.25,0.875,0.75,0.375,0.75,0.125,0.75,0.125");  // 001, binary metal-oxide prototype (ICSD #62764)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC32_40_c_3c"){
        vparameters.push_back("6.8124,1.7712407,0.89328284,0.57380348,0.61698844,0.21913764,0.85360935,0.57437175,0.93647361,0.96338832,0.74668906,0.23766452,0.14359267,0.53362863,0.35672424");  // 001, binary metal-oxide prototype (ICSD #16031)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oF40_43_b_ab"){
        vparameters.push_back("12.827,0.82014501,0.2992126,0.8107,0.4542,0.6329,0.5,0.0825,0.7334,0.1232");  // 001, binary metal-oxide prototype (ICSD #8014) //DX20210427 - equivalent to Ag2O3 (http://aflow.org/prototype-encyclopedia/A2B3_oF40_43_b_ab.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oI8_44_a_ac"){
        vparameters.push_back("4.931,0.61914419,1.0539444,0.4861,0.9779,0.2303,0.1279");  // 001, binary metal-oxide prototype (ICSD #180565)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP2_47_a_e"){
        vparameters.push_back("3.656,1.0765864,1.6829869");  // 001, binary metal-oxide prototype (ICSD #95729)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP3_47_a_m"){
        vparameters.push_back("1.3761,13.1138,4.7018385,0.67078142");  // 001, binary metal-oxide prototype (ICSD #54126)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_oP5_47_e_bdr"){
        vparameters.push_back("3.21,1.0183801,2.2482866,0.725");  // 001, binary metal-oxide prototype (ICSD #76022)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2_oP7_47_behq_ad"){
        vparameters.push_back("3.677,1.0598314,1.6907805,0.681");  // 001, binary metal-oxide prototype (ICSD #95462)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oP14_55_gh_bg"){
        vparameters.push_back("6.461,1.4519734,0.51335707,0.3108,0.0841,0.9497,0.17564,0.3014,0.8132");  // 001, binary metal-oxide prototype (ICSD #97282)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP12_59_e_2e"){
        vparameters.push_back("2.859,3.2402938,0.78926198,0.861,0.058,0.209,0.43,0.055,0.446");  // 001, binary metal-oxide prototype (ICSD #54114)
        vparameters.push_back("3.307,1.9525854,0.83979438,0.8896,0.4922,0.5744,0.7182,0.3388,0.951");  // 002, binary metal-oxide prototype (ICSD #83863)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP12_61_a_c"){
        vparameters.push_back("4.8,1.2666667,1.2520833,0.405,0.075,0.062");  // 001, binary metal-oxide prototype (ICSD #24774)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oP12_62_c_2c"){
        vparameters.push_back("9.2734,0.30881877,0.48762051,0.1332,0.9742,0.9655,0.2162,0.2796,0.6799");  // 001, binary metal-oxide prototype (ICSD #171866)
        vparameters.push_back("10.92,0.28296703,0.45549451,0.375,0,0.0417,0.25,0.2083,0.75");  // 002, binary metal-carbide prototype (ICSD #31973)
        vparameters.push_back("5.169,0.43180499,1.1501257,0.796,0.176,0.166,0.049,0.499,0.839");  // 003, binary metal-carbide prototype (ICSD #187138)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_oC20_63_ac_cf"){
        vparameters.push_back("2.676,3.3415546,2.6307922,0.246,0.8928,0.6343,0.0778");  // 001, binary metal-oxide prototype (ICSD #161062)
        vparameters.push_back("2.525,3.2772277,2.4554455,0.25,0.907,0.352,0.57");  // 002, binary metal-oxide prototype (ICSD #162252)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_oC28_63_cf_acf"){
        vparameters.push_back("2.6487,3.518443,3.538113,0.1127,0.45465,0.63397,0.92846,0.2835,0.88272");  // 001, binary metal-oxide prototype (ICSD #263010)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC6_65_a_j"){
        vparameters.push_back("3.8922,1.4786239,0.7816659,0.6232");  // 001, binary metal-oxide prototype (ICSD #180398)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oF12_69_a_h"){
        vparameters.push_back("12.4267,0.90617783,0.22034812,0.3847");  // 001, binary metal-oxide prototype (ICSD #85080)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oF12_69_a_g"){
        vparameters.push_back("12.382,0.90910192,0.22097399,0.386");  // 001, binary metal-oxide prototype (ICSD #150886)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI8_71_g_e"){
        vparameters.push_back("7.075,0.84565371,0.59378092,0.75,0.374");  // 001, binary metal-oxide prototype (ICSD #25528) //DX20210428 - equivalent to CsO (http://aflow.org/prototype-encyclopedia/AB_oI8_71_g_i.html, part 3)
        vparameters.push_back("5.44,0.88841912,0.671875,0.25,0.376");  // 002, binary metal-carbide prototype (ICSD #25705)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI32_72_2j_2j"){
        vparameters.push_back("14.5192,0.42786104,0.28964406,0.42918,0.83459,0.81013,0.16385,0.31001,0.66412,0.92918,0.33025");  // 001, binary metal-oxide prototype (ICSD #181462)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_oI28_74_ace_hi"){
        vparameters.push_back("5.945,0.99444912,1.4109336,0.625,0.9945,0.2346,0.2263,0.9945");  // 001, binary metal-oxide prototype (ICSD #31156)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI24_87_h_2h"){
        vparameters.push_back("9.815,0.29006623,0.83,0.65,0.795,0.84,0.542,0.84");  // 001, binary metal-oxide prototype (ICSD #20227)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI32_88_cd_f"){
        vparameters.push_back("6.883,1.3252942,0.763,0.3298,0.0372");  // 001, binary metal-oxide prototype (ICSD #202055)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP16_113_ef_e"){
        vparameters.push_back("7.39,0.52503383,0.2633,0.493,0.2569,0.9302,0.0021,0.7915,0.0241");  // 001, binary metal-oxide prototype (ICSD #86144)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP5_115_g_abc"){
        vparameters.push_back("3.92,1.4923469,0.223");  // 001, binary metal-oxide prototype (ICSD #168808)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP8_127_bg_a"){
        vparameters.push_back("5.2968,0.70636988,0.2368");  // 001, binary metal-oxide prototype (ICSD #77680)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tP4_129_c_a"){
        vparameters.push_back("4.3977,0.72685722,0.4");  // 001, binary metal-oxide prototype (ICSD #15301)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP6_129_c_ac"){
        vparameters.push_back("3.94,1.715736,0.182,0.672");  // 001, binary metal-oxide prototype (ICSD #6318)
        vparameters.push_back("3.96,1.6742424,0.7,0.27");  // 002, binary metal-oxide prototype (ICSD #631464)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP8_129_cd_c"){
        vparameters.push_back("5.2722,0.74356056,0.56,0.06");  // 001, binary metal-oxide prototype (ICSD #27961)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_tP16_130_cf_c"){
        vparameters.push_back("5.2759,1.4871775,0.00284,0.28474,0.02872");  // 001, binary metal-oxide prototype (ICSD #50732) //DX20210428 - equivalent to alpha-WO3 (http://aflow.org/prototype-encyclopedia/A3B_tP16_130_cf_c.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP10_132_i_be"){
        vparameters.push_back("5.84,0.97945205,0.236");  // 001, binary metal-oxide prototype (ICSD #168807)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI4_139_b_a"){
        vparameters.push_back("2.982,1.8051643");  // 001, binary metal-oxide prototype (ICSD #41617)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI32_140_ab_hk"){
        vparameters.push_back("8.597,0.8235431,0.213,0.725,0.925");  // 001, binary metal-oxide prototype (ICSD #28347)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tI24_141_h_c"){
        vparameters.push_back("3.784,2.5105708,0,0.2089");  // 001, binary metal-oxide prototype (ICSD #93098)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_tI28_141_cd_ae"){
        vparameters.push_back("5.818,1.7005844,0.617");  // 001, binary metal-oxide prototype (ICSD #77675)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_hR7_148_f_a"){
        vparameters.push_back("5.5599342,2.7519371,0.086102407,0.41928555,0.75270204");  // 001, binary metal-oxide prototype (ICSD #174039)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP5_150_d_e"){
        vparameters.push_back("3.93,1.5572519,0.25,0.75003333");  // 001, binary metal-oxide prototype (ICSD #26864)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP6_152_b_a"){
        vparameters.push_back("3.577,2.426894,0.54,0.255");  // 001, binary metal-oxide prototype (ICSD #24062)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR16_155_2c_def"){
        vparameters.push_back("5.5628914,5.6059222,0.37495,0.20823333,0.33332222,0.83325555,0.74993194,0.41662916,0.083388885");  // 001, binary metal-oxide prototype (ICSD #23402)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR2_160_a_a"){
        vparameters.push_back("3.8800376,2.5025788,0,0.275");  // 001, binary metal-oxide prototype (ICSD #30361)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hR10_160_ab_2b"){
        vparameters.push_back("8.0881082,1.2634305,0,0.522,0.03,0.265,0.656,0.763,0.158");  // 001, binary metal-oxide prototype (ICSD #168810)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_hP7_162_a_k"){
        vparameters.push_back("5.13,0.92397661,0.333,0.25");  // 001, binary metal-oxide prototype (ICSD #20042)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP8_162_k_d"){
        vparameters.push_back("5.318,0.93098909,0.31083333,0.2228");  // 001, binary metal-oxide prototype (ICSD #26557)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP8_162_ab_k"){
        vparameters.push_back("5.06,0.94466403,0.3333,0.236");  // 001, binary metal-oxide prototype (ICSD #23575)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_hP14_163_c_i"){
        vparameters.push_back("5.14,1.844358,0.33333333,0.99996667,0.375");  // 001, binary metal-oxide prototype (ICSD #17009)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP16_163_ac_i"){
        vparameters.push_back("5.15,1.8563107,0.33333333,0.33336667,0.118");  // 001, binary metal-oxide prototype (ICSD #20041)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_164_d_abd"){
        vparameters.push_back("3.214,1.2523335,0.7882,0.3134");  // 001, binary metal-oxide prototype (ICSD #189287)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_hR14_166_acd_ch"){
        vparameters.push_back("5.9276448,2.457949,0.625,0.75491,0.75492,0.23533");  // 001, binary metal-oxide prototype (ICSD #92356)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hR16_166_ab3c_4c"){
        vparameters.push_back("2.994405,19.999779,0.62439,0.20773,0.08456,0.85386,0.56158,0.68711,0.73077");  // 001, binary metal-oxide prototype (ICSD #181459)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hR8_167_b_e"){
        vparameters.push_back("5.6295232,2.7697869,0.59");  // 001, binary metal-oxide prototype (ICSD #27023)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_hP6_186_ab_b"){
        vparameters.push_back("3.1,2.683871,0.625,0.375,0");  // 001, binary metal-oxide prototype (ICSD #24923)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP12_191_gl_f"){
        vparameters.push_back("7.298,0.53425596,0.212");  // 001, binary metal-oxide prototype (ICSD #32001) //DX20210428 - equivalent to Hexagonal WO3 (http://aflow.org/prototype-encyclopedia/A3B_hP12_191_gl_f.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP8_194_f_f"){
        vparameters.push_back("3.1523,2.4461504,0.592,0.333");  // 001, binary metal-oxide prototype (ICSD #24143)
        vparameters.push_back("3.1525,1.222839,0.408,0.667");  // 002, binary metal-oxide prototype (ICSD #50658)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_cF56_216_abe_2e"){
        vparameters.push_back("8.3941,0.6253,0.8722,0.3817");  // 001, binary metal-oxide prototype (ICSD #65338)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cI24_217_c_abc"){
        vparameters.push_back("6.68,0.375,0.25");  // 001, binary metal-oxide prototype (ICSD #28387)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cP12_223_c_d"){
        vparameters.push_back("5.633");  // binary metal-oxide prototype (ICSD #181465)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_cP14_223_e_c"){
        vparameters.push_back("5.585");  // binary metal-oxide prototype (ICSD #30444)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_cP10_224_b_d"){
        vparameters.push_back("4.904");  // binary metal-oxide prototype (ICSD #15999) //DX20210428 - equivalent to Mg3P2 (D5_{5}, http://aflow.org/prototype-encyclopedia/A3B2_cP10_224_d_b.html, part 3)
      }
      // ---------------------------------------------------------------------------
      // ternaries 
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_aP5_1_a_3a_a"){
        vparameters.push_back("4.9154,1.0100094,1.1190951,118.0943,115.8701,60.69,0.3612,0.3496,0.571,0.9717,0.4206,0.7752,0.5628,0.2252,0.2364,0.2073,0.9587,0.2321,0.852,0.8545,0.0632");  // 001, ternary metal-oxide prototype (ICSD #9414)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_aP10_1_2a_2a_6a"){
        vparameters.push_back("4.8631,1.1422549,1.1471078,118.8837,87.8984,88.1302,0.0069620739,0.83060279,0.17002557,0.99323899,0.16946141,0.83022036,0.99843919,0.5017688,0.50237438,0.49892041,0.49733084,0.49721955,0.76260629,0.22750203,0.49446795,0.23810792,0.49941492,0.76762003,0.23680148,0.77225219,0.50516969,0.76138702,0.50020555,0.23221641,0.74799829,0.75478863,0.75712733,0.25364566,0.24534085,0.24238167");  // 001, ternary metal-oxide prototype (ICSD #92339)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_aP14_1_4a_8a_2a"){
        vparameters.push_back("5.07,1.0389349,1.6682446,90,90,118.7679,0.7051,0.3453,0.9454,0.3618,0.6583,0.0537,0.3609,0.6575,0.4456,0.7055,0.344,0.5538,0.3785,0.3817,0.9432,0.3782,0.3804,0.5554,0.6419,0.0012,0,0.6414,0,0.4988,0.9969,0.6188,0.0561,0.9962,0.6108,0.4425,0.467,0.6295,0.2505,0.8377,0.369,0.7502,0.0302,0.0599,0.75,0.9704,0.9402,0.25");  // 001, ternary metal-oxide prototype (ICSD #33532)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_aP14_1_4a_2a_8a"){
        vparameters.push_back("5.9816,0.99653939,0.99476729,60.6576,60.3581,60.4468,0.50008868,0.49963654,0.99969723,0.00080812292,0.00022962752,0.9993763,0.0012147868,0.4994824,0.99940545,0.0010131602,0.49938612,0.49974352,0.62714457,0.12088612,0.62306824,0.3779998,0.87504041,0.37413834,0.75869598,0.71308292,0.76798147,0.21826694,0.26068039,0.76005652,0.75773454,0.26821078,0.21444009,0.76002712,0.26116151,0.76095,0.78864778,0.73799057,0.23794214,0.23609507,0.28266831,0.24065217,0.23622713,0.74099405,0.78262508,0.23503613,0.7415005,0.24087342");  // 001, ternary metal-oxide prototype (ICSD #280061)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C3_aP16_1_2a_8a_6a"){
        vparameters.push_back("5.758,1.0804793,1.41438,111.2299,90,117.5651,0,0,0,0.6983,0,0.5002,0.2308,0.3415,0.4001,0.6044,0.1767,0.2096,0.2649,0.1746,0.7094,0.807,0.3421,0.8985,0.775,0.5403,0.6218,0.052,0.7323,0.3018,0.4597,0.5431,0.1192,0.3722,0.7319,0.8023,0.0141,0.3307,0.1411,0.0082,0.3312,0.6386,0.4408,0.0886,0.4538,0.3407,0.0912,0.9526,0.6739,0.6781,0.4173,0.6961,0.68,0.9185");  // 001, ternary metal-oxide prototype (ICSD #10473)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2_aP9_2_a_3i_i"){
        vparameters.push_back("9.168,0.38645288,0.70658813,87.75,110.34,88.12,0.9696,0.0027,0.2761,0.6574,0.0482,0.1104,0.6933,0.9972,0.5684,0.80721,0.01267,0.34537");  // 001, ternary metal-oxide prototype (ICSD #28151)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_aP10_2_i_3i_i"){
        vparameters.push_back("5.3998,1,1.0630023,62.948,62.948,60.4941,0.36208198,0.36208198,0.91375405,0.22928467,0.95503624,0.26123929,0.95503624,0.5544398,0.26123929,0.5544398,0.22928467,0.26123929,0.15135687,0.15135687,0.54592939");  // 001, ternary metal-oxide prototype (ICSD #29203)
        vparameters.push_back("3.592,1.3229399,2.3908686,88.32,79.6,89.3,0.28495,0.18523,0.57431,0.098,0.4968,0.3303,0.1939,0.0511,0.1509,0.7138,0.483,0.0906,0.1971,0.3855,0.13572");  // 002, ternary metal-oxide prototype (ICSD #82242)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C6_aP10_2_a_bi_3i"){
        vparameters.push_back("3.6909,1.5254003,1.9337289,90.3804,102.5515,91.1364,0.38911267,0.4711918,0.23520992,0.12490443,0.7954708,0.78978038,0.48088771,0.74389833,0.47547872,0.15767786,0.30795082,0.88739393");  // 001, ternary metal-oxide prototype (ICSD #411500)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_aP12_2_i_4i_i"){
        vparameters.push_back("4.6964,1.2410996,1.037731,88.37,87.56,82.79,0.5038,0.3379,0.2362,0.7519,0.6444,0.4171,0.7757,0.1231,0.4324,0.2707,0.6141,0.0956,0.2155,0.0964,0.0512,0.9834,0.8255,0.267");  // 001, ternary metal-oxide prototype (ICSD #4189)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8C2_aP13_2_ai_4i_i"){
        vparameters.push_back("5.196,1.0306005,1.2519246,110.78,91.31,68.08,0.7141,0.78382,0.31047,0.1399,0.687,0.3408,0.69,0.3711,0.1652,0.7635,0.8071,0.0267,0.614,0.8314,0.607,0.6288,0.6462,0.78267");  // 001, ternary metal-oxide prototype (ICSD #27184)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C3_aP14_2_i_3i_3i"){
        vparameters.push_back("5.7345,1.0271689,1.1073851,115.462,89.279,94.767,0.641,0.1907,0.6885,0.9371,0.7683,0.5556,0.1444,0.2345,0.9763,0.5634,0.2585,0.2555,0.754,0.1008,0.9005,0.6925,0.923,0.3758,0.7444,0.47,0.6571");  // 001, ternary metal-oxide prototype (ICSD #99580)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2_aP14_2_4i_i_abcf"){
        vparameters.push_back("6.6489,0.83507046,0.92004693,71.173,115.213,82.805,0.22,0.864,0.395,0.997,0.622,0.655,0.655,0.667,0.841,0.227,0.977,0.865,0.6406,0.6365,0.2529");  // 001, ternary metal-oxide prototype (ICSD #59657)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_aP16_2_4i_i_3i"){
        vparameters.push_back("7.0652,1.0085348,1.4067542,80.735,71.323,64.041,0.02,0.2563,0.4066,0.2064,0.2051,0.0107,0.2677,0.6523,0.1467,0.4423,0.2117,0.62,0.6844,0.2054,0.2365,0.132,0.012,0.676,0.453,0.18,0.207,0.73,0.447,0.177");  // 001, ternary metal-oxide prototype (ICSD #423336)
        vparameters.push_back("8.665,1.0240046,0.44893249,93.3,89.9,105,0.7377,0.4221,0.568,0.5317,0.2032,0.2275,0.8098,0.1951,0.3329,0.6938,0.2744,0.4141,0.8276,0.798,0.2116,0.7728,0.5487,0.6919,0.3962,0.1479,0.273,0.9113,0.1344,0.2669");  // 002, ternary metal-carbo-nitride prototype (ICSD #77046)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mP16_4_2a_2a_4a"){
        vparameters.push_back("5.415,1.9774515,0.58173592,90.886,0.807,0.114,0.464,0.693,0.386,0.464,0.275,0.13,0.911,0.225,0.37,0.911,0.026,0,0,0.526,0,0,0,0.251,0.958,0.5,0.251,0.958");  // 001, ternary metal-oxide prototype (ICSD #84868)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2_mC18_5_a_3c_c"){
        vparameters.push_back("9.18,0.38997821,0.70588235,79.6,0,0.03,0.002,0.278,0.336,0.003,0.117,0.308,0.019,0.551,0.1923,0.967,0.3477");  // 001, ternary metal-oxide prototype (ICSD #21067)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_mC18_5_3c_c_a"){
        vparameters.push_back("9.242,0.38151915,0.99167929,138.1526,0,0.457,0.949,0.109,0.314,0.975,0.284,0.376,0.489,0.568,0.538,0.983,0.347");  // 001, ternary metal-oxide prototype (ICSD #26998)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mC24_5_c_c_4c"){
        vparameters.push_back("12.52,0.30591054,0.98154153,148.8255,0.61644,0.99232,0.81355,0.12737,0,0.23054,0.88667,0.0063,0.1461,0.43654,0.9913,0.799,0.80607,0.9747,0.8609,0.37736,0.0239,0.5147");  // 001, ternary metal-oxide prototype (ICSD #429)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B3C2_mC26_5_4c_ac_c"){
        vparameters.push_back("7.46,0.82989276,1.2530831,116.63,0,0.248,0.809,0.943,0.23,0.203,0.982,0.693,0.951,0.709,0.459,0.984,0.763,0.113,0.051,0.719,0.29,0.524,0.591");  // 001, ternary metal-oxide prototype (ICSD #29359)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_mC28_5_c_4c_abc"){
        vparameters.push_back("10.565,0.60745859,0.76090866,88.95,0,0.9494,0.6623,0.9743,0.6919,0.6911,0.051,0.4901,0.5627,0.767,0.7127,0.8026,0.92,0.7833,0.5837,0.1848,0.7993,0.33358,0.9929,0.82233");  // 001, ternary metal-oxide prototype (ICSD #280608)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mP5_6_b_2ab_a"){
        vparameters.push_back("3.9729,0.99723124,0.99723124,90,0.0052,0.4792,0.5208,0.9948,0.4907,0.5093,0,0,0.5085,0.4915");  // 001, ternary metal-oxide prototype (ICSD #186460)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_mC12_8_2a_a_3a"){
        vparameters.push_back("4.8339,0.94296531,2.3231138,107.0815,0.415,0.727,0.099,0.281,0.803,0.505,0.314,0.497,0.62,0.338,0.84,0.677");  // 001, ternary metal-oxide prototype (ICSD #95439)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C5_mC22_8_a_5a_5a"){
        vparameters.push_back("9.964,0.41629868,0.97859293,144.9445,0,0,0.554,0.543,0.845,0.138,0.182,0.852,0.376,0.663,0.628,0.345,0.57,0.071,0.182,0.336,0.86,0.693,0.699,0.838,0.342,0.191");  // 001, ternary metal-oxide prototype (ICSD #203031)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC24_8_ab_ab_2a2b"){
        vparameters.push_back("5.005,1.7322677,1.0021978,109.428,0,0,0.5,0.5,0.746604,0.239812,0.227089,0.733288,0,0.3333,0,0.5,0.3333,0.5,0.760083,0.3333,0.280249,0.753099,0.157997,0.733288");  // 001, ternary metal-oxide prototype (ICSD #164853)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C5_mC26_8_2a_2a2b_a2b"){
        vparameters.push_back("6.8558,1.6648677,0.96866886,86.842,0.8454,0.1356,0.4542,0.5013,0.0033,0.6425,0.2951,0.9849,0.5812,0.2367,0.6476,0.65634,0.3166,0.647,0.78703,0.8162,0.35,0.1383,0.593,0.444,0.3599,0.044");  // 001, ternary metal-oxide prototype (ICSD #174312)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C5_mC26_8_2a2b_2a_a2b"){
        vparameters.push_back("7.4722,1.6430235,0.97827949,83.34,0,0.3671,0.705,0,0.1591,0.8645,0.5413,0.4999,0.396,0.726,0.3523,0.75462,0.6824,0.3548,0.78755,0.1797,0.638,0.126,0.42,0.556,0.368,0.962");  // 001, ternary metal-oxide prototype (ICSD #174313)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mC20_9_a_a_3a"){
        vparameters.push_back("9.6664,0.57731937,0.82061574,144.3355,0,0,0,0.44042,0,0.66063,0.1892,0.2675,0.8554,0.12705,0.79465,0.8554,0.38915,0.93785,0.8446");  // 001, ternary metal-oxide prototype (ICSD #247765)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_mC24_9_2a_3a_a"){
        vparameters.push_back("5.427,1.6629814,1.107518,123.625,0.5,0.83,0,0.5,0.5,0,0.486,0.989,0.77,0.549,0.695,0.784,0.966,0.842,0.75,0.5,0.1584,0");  // 001, ternary metal-oxide prototype (ICSD #31941)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C3_mC32_9_a_4a_3a"){
        vparameters.push_back("10.967,0.53068296,0.97586396,135.4126,0.2656,0.1036,0.2441,0.4238,0.0892,0.631,0.5448,0.4621,0.351,0.3323,0.588,0.445,0.6464,0.909,0.532,0.2931,0.345,0.109,0.6962,0.261,0.843,0.3569,0.806,0.296");  // 001, ternary metal-oxide prototype (ICSD #1410)
        vparameters.push_back("10.993,0.52296916,0.97277358,135.9926,0.7199,0.099,0.7082,0.625,0.096,0.342,0.447,0.462,0.6,0.656,0.593,0.509,0.352,0.929,0.429,0.702,0.247,0.87,0.804,0.234,0.641,0.622,0.803,0.654");  // 002, ternary metal-oxide prototype (ICSD #10501)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C_mC32_9_4a_3a_a"){
        vparameters.push_back("5.8277,2.8612317,1.1510373,124.4407,0.2762,0.5042,0.7585,0.2333,0.7169,0.7495,0.2426,0.2932,0.7216,0.6423,0.3917,0.6905,0.7546,0.3167,0.0123,0.7177,0.3137,0.4545,0.1733,0.4115,0.4685,0.7615,0.3933,0.2521");  // 001, ternary metal-oxide prototype (ICSD #49624)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mP10_11_a_e_d2e"){
        vparameters.push_back("3.9754,1.9854606,0.99069276,89.223,0.4833484,0.5081843,0.43760642,0.012023374,0.97729189,0.5328779");  // 001, ternary metal-oxide prototype (ICSD #162834)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mP12_11_e_aef_e"){
        vparameters.push_back("7.2255,0.80074735,0.55788527,106,0.7293,0.348,0.1721,0.5932,0.2003,0.1482,0.3769,0.4781,0.1944");  // 001, ternary metal-oxide prototype (ICSD #155516)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_mP14_11_2e_e_4e"){
        vparameters.push_back("7.056,0.53326247,0.78741497,101.36,0.871,0.8454,0.8882,0.3434,0.4273,0.7309,0.6934,0.0612,0.6968,0.5472,0.0746,0.6377,0.0713,0.1361");  // 001, ternary metal-oxide prototype (ICSD #172780)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC8_12_a_c_i"){
        vparameters.push_back("6.3541,0.43247667,0.96973608,121.6413,0.425,0.77");  // 001, ternary metal-oxide prototype (ICSD #15098)
        vparameters.push_back("5.63,0.5079929,1.1192718,122.4891,0.483,0.205");  // 002, ternary metal-oxide prototype (ICSD #16270)
        vparameters.push_back("5.733,0.47402756,0.98044654,120.6617,0.5563,0.2433");  // 003, ternary metal-oxide prototype (ICSD #74978)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC8_12_c_a_i"){
        vparameters.push_back("5.33,0.53658537,1.1684803,122.7843,0.483,0.205");  // 001, ternary metal-oxide prototype (ICSD #26609)
        vparameters.push_back("5.53,0.52151899,1.2648463,125.3142,0.727,0.143");  // 002, ternary metal-oxide prototype (ICSD #30379)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC8_12_a_d_i"){
        vparameters.push_back("5.933,0.58183044,0.99039272,88.74,0.6657,0.7933");  // 001, ternary metal-oxide prototype (ICSD #74848)
        vparameters.push_back("6.0756,0.46230825,0.96662058,72.013,0.1457,0.301");  // 002, ternary metal-oxide prototype (ICSD #95089)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_mC8_12_a_i_b"){
        vparameters.push_back("6.7317,0.69509931,0.99175543,155.5107,0.081292291,0.39561177");  // 001, ternary metal-oxide prototype (ICSD #420136)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_mC8_12_a_i_c"){
        vparameters.push_back("5.1758,0.58112756,1.2351327,122.8692,0.5379,0.3011");  // 001, ternary metal-oxide prototype (ICSD #420138)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_mC10_12_a_i_i"){
        vparameters.push_back("5.6346,0.4943208,0.87376211,82.872,0.3079,0.1182,0.3684,0.7583");  // 001, ternary metal-oxide prototype (ICSD #174134)
        vparameters.push_back("5.7877,0.98536552,1.1609793,125.5059,0.6126,0.2855,0.937,0.7903");  // 002, ternary metal-carbo-nitride prototype (ICSD #411094)
        vparameters.push_back("5.0475,0.99078752,1.2039624,121.1857,0.07199,0.22779,0.36559,0.69161");  // 003, ternary metal-carbo-nitride prototype (ICSD #411341)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC16_12_ad_i_2i"){
        vparameters.push_back("10.543,0.27151665,0.47091909,80.6455,0.2491,0.221,0.381,0.867,0.14,0.615");  // 001, ternary metal-oxide prototype (ICSD #165326)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C4_mC18_12_a_2i_j"){
        vparameters.push_back("9.1727,0.68986231,0.66786224,74.5,0.8075,0.5892,0.6212,0.1423,0.6086,0.7113,0.799");  // 001, ternary metal-oxide prototype (ICSD #47223)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C4_mC18_12_i_ai_2i"){
        vparameters.push_back("9.946,0.27930826,0.90709833,135.3216,0.612,0.7729,0.27,0.632,0.2771,0.8442,0.0647,0.3211");  // 001, ternary metal-oxide prototype (ICSD #66509)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C2_mC18_12_i_a2i_i"){
        vparameters.push_back("11.37,0.33421284,0.58223395,79.9,0.4022,0.3561,0.128,0.336,0.325,0.992,0.1495,0.0928");  // 001, ternary metal-oxide prototype (ICSD #36097) //DX20210427 - equivalent to K2O5Ti2 (http://aflow.org/prototype-encyclopedia/A2B5C2_mC18_12_i_a2i_i.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4_mC18_12_2i_a_j"){
        vparameters.push_back("10.3166,0.65099936,0.64017215,71.882,0.8093,0.5945,0.6152,0.1611,0.6051,0.7017,0.8181");  // 001, ternary metal-oxide prototype (ICSD #72287)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mC20_12_i_ghi_e"){
        vparameters.push_back("8.5311,0.99410393,0.71146745,135.1955,0.785,0.7235,0.249,0.4958,0.745,0.0524");  // 001, ternary metal-oxide prototype (ICSD #94313)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_mC22_12_b_ac3i_i"){
        vparameters.push_back("20.1163,0.18255345,1.0184378,168.6098,0.14373409,0.45200889,0.58162295,0.98359569,0.79223864,0.98296547,0.53586962,0.85028618");  // 001, ternary metal-oxide prototype (ICSD #28471)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C_mC22_12_i_2ij_a"){
        vparameters.push_back("9.8821,0.60411249,0.5453092,87.182,0.3297,0.691,0.1281,0.761,0.3793,0.418,0.6136,0.756,0.143");  // 001, ternary metal-oxide prototype (ICSD #280433)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C2_mC22_12_aij_i_h"){
        vparameters.push_back("6.9324,1.217385,0.72595349,71.728,0.6842,0.5909,0.205,0.7826,0.9049,0.7849,0.8481,0.7085");  // 001, ternary metal-oxide prototype (ICSD #250002)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC8_mC24_12_ci_a_2ij"){
        vparameters.push_back("8.593,0.63609915,0.88921215,83.35,0.625,0.707,0.225,0.449,0.682,0.892,0.524,0.25,0.685");  // 001, ternary metal-oxide prototype (ICSD #15898)
        vparameters.push_back("8.4945,0.64468774,0.80041203,88.588,0.636,0.7291,0.2256,0.4426,0.7118,0.9384,0.5235,0.2539,0.7139");  // 002, ternary metal-oxide prototype (ICSD #82620)
        vparameters.push_back("8.4711,0.66584033,0.77112772,90,0.3576,0.2641,0.7771,0.5478,0.7428,0.9571,0.9721,0.7614,0.2741");  // 003, ternary metal-oxide prototype (ICSD #155504)
        vparameters.push_back("8.8965,0.61409543,0.88498848,86.507,0.3615,0.2958,0.7909,0.5699,0.316,0.0984,0.4638,0.2442,0.3289");  // 004, ternary metal-oxide prototype (ICSD #246524)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mC24_12_i_i_4i"){
        vparameters.push_back("12.308,0.30638609,0.9808661,148.7576,0.6041,0.7988,0.1268,0.2291,0.884,0.146,0.437,0.799,0.811,0.869,0.376,0.511");  // 001, ternary metal-oxide prototype (ICSD #14016) //DX20210428 - equivalent to AlNbO4 (http://aflow.org/prototype-encyclopedia/ABC4_mC24_12_i_i_4i.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mC24_12_i_2ij_g"){
        vparameters.push_back("9.069,0.62906605,0.50071673,87.71,0.2284,0.741,0.4858,0.8859,0.7791,0.5977,0.1975,0.3496,0.2551,0.3012");  // 001, ternary metal-oxide prototype (ICSD #4164)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_mC24_12_acg_h_ij"){
        vparameters.push_back("4.8984,1.7249714,1.1695656,124.9135,0.3132,0.3296,0.4728,0.7125,0.4799,0.6679,0.7288");  // 001, ternary metal-oxide prototype (ICSD #153094)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C_mC24_12_2i_b3i_a"){
        vparameters.push_back("13.469,0.29892345,0.72647561,135.4662,0.593,0.846,0.195,0.6111,0.0572,0.6899,0.2001,0.1126,0.3499,0.5538");  // 001, ternary metal-oxide prototype (ICSD #65032)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B6C_mC24_12_dgh_ij_a"){
        vparameters.push_back("5.0679,1.722903,0.99238343,69.76,0.3314,0.1762,0.7756,0.7723,0.7706,0.1522,0.2361");  // 001, ternary metal-oxide prototype (ICSD #38381)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mC24_12_i_2ij_g"){
        vparameters.push_back("9.3483,0.60942631,0.24667587,88.01,0.2263,0.747,0.01,0.8751,0.564,0.5896,0.49,0.3517,0.2697,0.5866");  // 002, ternary metal-oxide prototype (ICSD #173875)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B6C_mC24_12_dgh_ij_a"){
        vparameters.push_back("5.6877,1.5351372,0.98588182,68.99,0.3333,0.1737,0.7974,0.7994,0.7929,0.1377,0.2129");  // 002, ternary metal-oxide prototype (ICSD #38382)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8C_mC24_12_ci_2ij_a"){
        vparameters.push_back("8.9606,0.61453474,0.86301141,87.16,0.6396,0.2942,0.2082,0.5729,0.3166,0.9082,0.537,0.7459,0.3271");  // 001, ternary metal-oxide prototype (ICSD #155505)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C8_mC24_12_a_ci_2ij"){
        vparameters.push_back("8.641,0.63894225,0.83270455,86.2,0.3444,0.7152,0.784,0.4475,0.7473,0.0781,0.4566,0.7447,0.6936");  // 001, ternary metal-oxide prototype (ICSD #155507)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC8_mC24_12_di_a_2ij"){
        vparameters.push_back("9.8811,0.57525984,0.83231624,141.7397,0.2458,0.6199,0.4712,0.7135,0.9624,0.6953,0.7335,0.7542,0.7248");  // 001, ternary metal-oxide prototype (ICSD #155851)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC8_mC26_12_2i_a_4i"){
        vparameters.push_back("13.8199,0.21280907,0.70710352,135,0.6648,0.5127,0.8481,0.1833,0.7041,0.3592,0.6524,0.9503,0.6629,0.7054,0.0419,0.3801");  // 001, ternary metal-oxide prototype (ICSD #491)
        vparameters.push_back("10.0206,0.31462188,1.4165319,134.9061,0.476,0.8274,0.8068,0.6379,0.666,0.807,0.067,0.849,0.263,0.812,0.604,0.461");  // 002, ternary metal-oxide prototype (ICSD #80336)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C8_mC26_12_i_ah_2ij"){
        vparameters.push_back("9.695,0.58122744,0.50665291,76.65,0.7411,0.2803,0.0693,0.3964,0.3887,0.1021,0.6366,0.6188,0.2728,0.1335");  // 001, ternary metal-oxide prototype (ICSD #971)
        vparameters.push_back("10.806,0.53747918,0.95049972,153.0887,0.239,0.5632,0.843,0.683,0.601,0.167,0.566,0.201,0.767,0.11");  // 002, ternary metal-oxide prototype (ICSD #16957)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8C2_mC26_12_ai_2ij_i"){
        vparameters.push_back("9.6715,0.56173293,0.93373313,135.7332,0.3995,0.5998,0.0268,0.2932,0.6325,0.9584,0.8106,0.218,0.7119,0.2314,0.6942");  // 001, ternary metal-oxide prototype (ICSD #412273)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C5_mC26_12_ghi_i_aj"){
        vparameters.push_back("7.474,1.6440995,0.97792347,83.59,0.7547,0.2125,0.3522,0.3161,0.1918,0.8175,0.7084,0.1295,0.2715");  // 001, ternary metal-oxide prototype (ICSD #73134)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C6_mC26_12_2i_ai_hj"){
        vparameters.push_back("8.9953,0.67839872,0.86572988,88.335,0.8,0.32054,0.37493,0.3507,0.89019,0.03725,0.66914,0.57061,0.2988,0.8419");  // 001, ternary metal-oxide prototype (ICSD #426450)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C2_mC26_12_aci_bij_i"){
        vparameters.push_back("10.222,0.60955782,0.71238505,78.69,0.675,0.2878,0.2219,0.1443,0.3672,0.1986,0.3865,0.2156,0.3265");  // 001, ternary metal-oxide prototype (ICSD #250388)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C8_mC26_12_a_2i_4i"){
        vparameters.push_back("10.006,0.28642814,0.97401559,88.83,0.33595,0.34692,0.14815,0.8317,0.2944,0.1495,0.342,0.8019,0.3452,0.5431,0.9582,0.8231");  // 001, ternary metal-oxide prototype (ICSD #62096)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B8C_mC26_12_2i_4i_a"){
        vparameters.push_back("11.8,0.26677966,0.84032203,129.5074,0.6049,0.3627,0.1694,0.7145,0.9492,0.7829,0.8187,0.4652,0.3646,0.9062,0.4168,0.3112");  // 001, ternary metal-oxide prototype (ICSD #155638)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C3_mC28_12_i_ghi_ij"){
        vparameters.push_back("7.09,1.5669958,0.91396333,78.4,0.7379,0.8106,0.8318,0.144,0.6265,0.6468,0.0999,0.1628,0.7079,0.8634,0.2625");  // 001, ternary metal-oxide prototype (ICSD #2269)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_mC28_12_i_2i_4i"){
        vparameters.push_back("13.552,0.26630756,0.86739965,133.88,0.5485,0.163,0.6581,0.5391,0.8591,0.1989,0.745,0.278,0.647,0.945,0.882,0.616,0.429,0.349");  // 001, ternary metal-oxide prototype (ICSD #69734)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_mC28_12_2i_i_2ij"){
        vparameters.push_back("12.348,0.49246842,0.91784095,143.1937,0.7258,0.2352,0.4097,0.255,0.05113,0.22671,0.0491,0.3824,0.274,0.3618,0.4355,0.2666,0.0785");  // 001, ternary metal-oxide prototype (ICSD #16154)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mC28_12_2i_4i_i"){
        vparameters.push_back("19.3,0.22580311,0.31635751,85.199,0.9025,0.8237,0.5785,0.7212,0.601,0.07,0.693,0.61,0.778,0.87,0.901,0.22,0.7608,0.253");  // 001, ternary metal-oxide prototype (ICSD #80668)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B10C2_mC30_12_ah_3ij_i"){
        vparameters.push_back("7.911,1.0065731,0.88054608,64.24,0.25681,0.6394,0.1059,0.09758,0.2814,0.3214,0.4769,0.36603,0.17838,0.36007,0.76265,0.21221");  // 001, ternary metal-oxide prototype (ICSD #10104)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B3C2_mC30_12_3ij_ah_i"){
        vparameters.push_back("10.985,0.51297223,0.58734638,74.7,0.7701,0.1918,0.884,0.3839,0.525,0.8611,0.541,0.67081,0.8545,0.4807,0.2542,0.8085");  // 001, ternary metal-oxide prototype (ICSD #50707)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_mC32_12_i_g2ij_2i"){
        vparameters.push_back("14.0076,0.40295982,0.90732174,149.66,0.746,0.4913,0.6615,0.6845,0.5934,0.4052,0.7013,0.0292,0.8854,0.2511,0.7308,0.628,0.7629,0.7757");  // 001, ternary metal-oxide prototype (ICSD #29269)
        vparameters.push_back("14.018,0.40540733,0.90843202,149.5084,0.75,0.4815,0.653,0.3489,0.4207,0.3459,0.6297,0.964,0.1075,0.7425,0.266,0.6157,0.7626,0.7809");  // 002, ternary metal-oxide prototype (ICSD #201562)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C_mC32_12_g2ij_2i_i"){
        vparameters.push_back("14.2188,0.40801615,0.92097083,149.1252,0.8269,0.6806,0.6033,0.4145,0.7362,0.019,0.8821,0.2543,0.7428,0.4891,0.6589,0.6194,0.7534,0.7805");  // 001, ternary metal-oxide prototype (ICSD #61399)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mP12_13_e_2g_f"){
        vparameters.push_back("6.7723,0.84885489,0.72418528,136.0686,0.19564,0.31412,0.23115,0.64564,0.9186,0.21435,0.06685,0.88984");  // 001, ternary metal-oxide prototype (ICSD #236416)
        vparameters.push_back("7.0885,0.78957466,0.7656768,132.7272,0.23456,0.31286,0.7432,0.0942,0.3529,0.2776,0.4376,0.7824");  // 002, ternary metal-oxide prototype (ICSD #415430)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_mP16_13_e_e2f_2g"){
        vparameters.push_back("5.871,1.1405212,0.96389031,109.8,0.864,0.381,0.126,0.606,0.797,0.103,0,0.768,0.342,0.529");  // 001, ternary metal-oxide prototype (ICSD #10319)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3_mP14_14_e_e_ae"){
        vparameters.push_back("6.892,0.97547882,0.9556007,123.1916,0.33,0.8757,0.2624,0.1809,0.1498,0.7585,0.3984,0.6144,0.1311");  // 001, ternary metal-oxide prototype (ICSD #16223)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_mP16_14_e_2e_e"){
        vparameters.push_back("6.277,0.98534332,0.98932611,121.3071,0.25,0.375,0.25,0,0.85,0.8,0.45,0.9,0.7,0.25,0.875,0.25");  // 001, ternary metal-oxide prototype (ICSD #45511)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mP16_14_e_e_2e"){
        vparameters.push_back("5.824,1.135989,0.9946772,102.48,0.7992,0.06802,0.6937,0.351,0.138,0.842,0.0859,0.2612,0.5114,0.5705,0.1275,0.2863");  // 001, ternary metal-oxide prototype (ICSD #47116)
        vparameters.push_back("6.8521,1.018622,0.81020417,113.212,0.768668,0.119738,0.565776,0.24219,0.13389,0.005089,0.78199,0.126728,0.815387,0.700386,0.122027,0.346386");  // 002, ternary metal-carbo-nitride prototype (ICSD #412278)
        vparameters.push_back("4.401,2.18057259714,1.58638945694,108.274,0.4105,0.1407,0.2101,-0.0732,0.4107,0.1839,0.2034,0.4127,0.3067,0.8385,0.2999,0.0865");  // 003, (part 3)
        vparameters.push_back("5.304,0.994909502262,1.0,114.38,0.284,0.027,0.725,0.76316,0.01033,0.75464,0.3749,0.1238,0.6279,0.8752,0.1256,0.1206");  // 004, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC16_15_e_e_f"){
        vparameters.push_back("7.394,0.98214769,0.812429,129.1818,0.6336,0.1152,0.298,0.679,0.755");  // 001, ternary metal-oxide prototype (ICSD #10317)
        vparameters.push_back("7.826,1.0092001,0.84209047,131.9936,0.6053,0.1334,0.302,0.934,0.274");  // 002, ternary metal-oxide prototype (ICSD #16239)
        vparameters.push_back("8.274,1.1124003,0.86337926,135.1074,0.83797,0.39059,0.2843,0.8031,0.7094");  // 003, ternary metal-oxide prototype (ICSD #406564)
        vparameters.push_back("7.084,1.595567476,0.761434217956,113.2,-0.00332,0.3572,0.196,0.1098,0.1068");  // 004, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_mC16_15_e_f_e"){
        vparameters.push_back("8.062,1.0405358,0.85297693,133.6306,0.65223,0.1171,0.2167,0.8097,0.286");  // 001, ternary metal-oxide prototype (ICSD #407208)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mC20_15_a_e_ef"){
        vparameters.push_back("9.15,0.57897268,0.81843716,144.6223,0.5,0.0304,0.2329,0.7329,0.7177");  // 001, ternary metal-oxide prototype (ICSD #180415)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_mC20_15_ac_f_e"){
        vparameters.push_back("8.8523,0.93641201,0.68535861,137.3424,0.3775,0.2107,0.4333,0.1315");  // 001, ternary metal-oxide prototype (ICSD #400657)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mC24_15_a_e_2f"){
        vparameters.push_back("11.282,0.53669562,0.92867399,152.9206,0.56199,0.8162,0.1594,0.468,0.3953,0.8859,0.0312");  // 001, ternary metal-oxide prototype (ICSD #2533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mC24_15_e_e_2f"){
        vparameters.push_back("5.5735,2.0708352,0.93583924,85.93,0.6294,0.1037,0.2359,0.9661,0.9479,0.1442,0.7964,0.5111");  // 001, ternary metal-oxide prototype (ICSD #10123)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mC24_15_e_2f_e"){
        vparameters.push_back("7.037,1.5553503,0.75287765,134.07,0.1445,0.6212,0.7442,0.5418,0.3374,0.7919,0.7819,0.2968");  // 001, ternary metal-oxide prototype (ICSD #20335)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mC24_15_a_2f_e"){
        vparameters.push_back("11.459,0.50885767,0.90945981,153.6578,0.5586,0.6876,0.6574,0.0328,0.1153,0.373,0.4845");  // 001, ternary metal-oxide prototype (ICSD #169668)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_mC28_15_f_e_2f"){
        vparameters.push_back("8.7352,1.322454,0.91543411,146.0878,0.85025,0.08681,0.64286,0.97811,0.2354,0.7008,0.3515,0.1906,0.0562,0.799");  // 001, ternary metal-oxide prototype (ICSD #90084)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_mC32_15_cf_2f_e"){
        vparameters.push_back("10.1885,0.48830544,1.0641508,122.0653,0.2138,0.10808,0.21018,0.59805,0.1571,0.9845,0.7886,0.5538,0.9159,0.3971");  // 001, ternary metal-oxide prototype (ICSD #249417)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_mC32_15_f_e2f_c"){
        vparameters.push_back("10.101,0.4986635,0.99368379,139.175,0.5484,0.3424,0.25,0.7273,0.2807,0.0933,0.8417,0.1115,0.0913,0.4456");  // 001, ternary metal-oxide prototype (ICSD #24416)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC32_15_2e_2e_2f"){
        vparameters.push_back("6.5992,1.9896351,0.91223179,124.3691,0.054,0.312,0.562,0.812,0.24,0.9375,0.23,0.25,0.688,0.2");  // 001, ternary metal-oxide prototype (ICSD #2739)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_oP16_31_ab_2ab_a"){
        vparameters.push_back("6.3259,0.86090517,0.7820073,0.8326,0.4848,0.8704,0.8952,0.1736,0.3478,0.17039,0,0.753,0.6685,0.9872,0.7761,0.3196,0.891");  // 001, ternary metal-oxide prototype (ICSD #19002)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC20_36_a_a_3a"){
        vparameters.push_back("3.1364,3.1511924,2.3271585,0.25049,0.01,0.00141,0.7618,0.92506,0.0083,0.62985,0.809,0.62443,0.2085");  // 001, ternary metal-oxide prototype (ICSD #159030)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oC20_36_a_ab_a"){
        vparameters.push_back("5.243,1.9126454,1.0089643,0.42673,0,0.233,0.975,0.0686,0.989,0.7617,0.995,0.736");  // 001, ternary metal-oxide prototype (ICSD #284)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oC20_36_2a_a_2a"){
        vparameters.push_back("2.82,3.5957447,2.9361702,0.1479,0.3306,0.3944,0.1526,0.1219,0,0.4799,0.4099,0.2795,0.5922");  // 001, ternary metal-oxide prototype (ICSD #14159)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_36_a_4a_a"){
        vparameters.push_back("3.9313,3.7628774,1.4290438,0.1699,0.843,0.304,0.592,0.331,0.129,0.47,0.468,0.915,0.75,0.4145,0.805");  // 001, ternary metal-oxide prototype (ICSD #97688)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_oC32_36_a_5a_2a"){
        vparameters.push_back("3.607,2.6897699,2.9564735,0.817,0.217,0.765,0.635,0.538,0.782,0.708,0.365,0.427,0.431,0.57,0.012,0.6222,0.498,0.9342,0.5977");  // 001, ternary metal-oxide prototype (ICSD #23817)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC10_38_a_b_ae"){
        vparameters.push_back("3.971,1.4333921,1.4401914,0.0138,0.5364,0.5,0.7476,0.7842");  // 001, ternary metal-oxide prototype (ICSD #9533)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oC10_38_b_ad_a"){
        vparameters.push_back("3.9594,1.4210739,1.4253422,0.491,0.985,0,0.736,0.752");  // 001, ternary metal-oxide prototype (ICSD #154346)
        vparameters.push_back("4.0094,1.4020552,1.4063451,0.99,0.51,0,0.2525,0.2396");  // 002, ternary metal-oxide prototype (ICSD #161341)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oC10_38_a_ae_b"){
        vparameters.push_back("4.0094,1.4020552,1.4063451,0,0.49,0.01,0.2525,0.2396");  // 001, ternary metal-oxide prototype (ICSD #161419)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_40_b_a_b"){
        vparameters.push_back("5.086,2.0129768,1.1598506,0.0062,0.8275,0.5,0.8671,0");  // 001, ternary metal-oxide prototype (ICSD #37079)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oC28_40_c_abc_b"){
        vparameters.push_back("8.458,1.0917475,0.59943249,0.3588,0.8151,0.3473,0.9701,0,0.0541,0.8283,0.4673,0.0564,0.8092,0.8123");  // 001, ternary metal-oxide prototype (ICSD #80128)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_oI20_46_a2b_a_b"){
        vparameters.push_back("7.8341,0.70251848,0.70260783,0.754,0.2492,0.7247,0.5285,0.2246,0.4778,0.50044,0.24776");  // 001, ternary metal-oxide prototype (ICSD #182248)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP4_47_a_g_s"){
        vparameters.push_back("3.005,1.1930116,1.9464226,0.753");  // 001, ternary metal-oxide prototype (ICSD #15115)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_oP5_47_bce_h_a"){
        vparameters.push_back("3.875,1.0322581,1.0867097");  // 001, ternary metal-oxide prototype (ICSD #27949)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP12_50_a_m_c"){
        vparameters.push_back("5.3413,1.0090053,1,0.46497053,0.535959,0.74433453");  // 001, ternary metal-oxide prototype (ICSD #164740)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2_oP14_55_gh_a_h"){
        vparameters.push_back("6.159,1.6363046,0.5685988,0.64,0.69,0.78,0.95,0.923,0.681");  // 001, ternary metal-oxide prototype (ICSD #4418)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oP14_55_h_gh_a"){
        vparameters.push_back("5.546,1.7829066,0.57573026,0.63,0.69,0.942,0.68,0.79,0.95");  // 001, ternary metal-oxide prototype (ICSD #9010)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP8_58_c_g_a"){
        vparameters.push_back("4.967,1.0191262,0.55788202,0.3365,0.2555");  // 001, ternary metal-oxide prototype (ICSD #48007)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_oP12_59_e_a_be"){
        vparameters.push_back("4,2.3375,0.7975,0.631,0.453,0.914,0.82,0.094,0.131");  // 001, ternary metal-oxide prototype (ICSD #4202)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oP12_59_e_be_a"){
        vparameters.push_back("3.8831,2.7352888,0.84442842,0.755,0.387,0.928,0.72,0.064,0.247");  // 001, ternary metal-oxide prototype (ICSD #95957)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_oP16_59_b_a2e_e"){
        vparameters.push_back("3.61,3.133241,1.3296399,0.4805,0.8592,0.11452,0.94197,0.92698,0.51231,0.09788,0.60781");  // 001, ternary metal-oxide prototype (ICSD #59345)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_59_e_2e_e"){
        vparameters.push_back("3.42,3.2982456,0.77631579,0.92,0.86,0.294,0.642,0.977,0.506,0.141,0.086");  // 001, ternary metal-oxide prototype (ICSD #21013)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_oP16_59_a_b2e_e"){
        vparameters.push_back("3.5732,3.1778798,1.3026979,0.28,0.9716,0.6221,0.466,0.424,0,0.5999,0.1152");  // 001, ternary metal-oxide prototype (ICSD #88643)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_59_e_e_2e"){
        vparameters.push_back("3.448,3.2888631,0.77233179,0.863,0.1,0.09,0.82,0.2,0.44,0.03,0.36");  // 001, ternary metal-oxide prototype (ICSD #27769)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oP16_62_c_c_2c"){
        vparameters.push_back("11.409,0.3045841,0.46752564,0.359739,0.95718,0.0666,0.9087,0.0276,0.2569,0.2024,0.7013");  // 001, ternary metal-oxide prototype (ICSD #422560)
        vparameters.push_back("9.4253,0.301868375542,0.466902910252,0.1445,-0.04472,0.80108,0.28766,-0.05338,0.80286,-0.0876,0.5905");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC12_63_c_a_c"){
        vparameters.push_back("5.899,1.7355484,0.86218003,0.3285,0.8739");  // 001, ternary metal-oxide prototype (ICSD #40161)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_c_f"){
        vparameters.push_back("4.37,2.6887872,1.2402746,0,0.7,0.106,0");  // 001, ternary metal-oxide prototype (ICSD #15095)
        vparameters.push_back("5.22,2.335249,1.0402299,0.681,0,0.103,0");  // 002, ternary metal-oxide prototype (ICSD #15097)
        vparameters.push_back("3.8956,2.9194476,1.4593387,0.016,0.688,0.118,0");  // 003, ternary metal-oxide prototype (ICSD #409547)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_a_g"){
        vparameters.push_back("5.735,1.6024412,0.82911944,0.655,0.25,0.92");  // 001, ternary metal-oxide prototype (ICSD #15760)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC20_63_a_c_cf"){
        vparameters.push_back("3.0397,3.4074415,2.3609567,0.7493,0.06234,0.60159,0.9432");  // 001, ternary metal-oxide prototype (ICSD #261371)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_63_c_acf_c"){
        vparameters.push_back("3.839,3.8848659,1.4415212,0.688,0.409,0.9161,0.188,0.515");  // 001, ternary metal-oxide prototype (ICSD #10481)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC24_63_g_ce_c"){
        vparameters.push_back("11.3074,0.69348391,0.53453491,0.6501,0.9115,0.8768,0.8369,0.3537");  // 001, ternary metal-oxide prototype (ICSD #62140)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oC24_63_c_a_fg"){
        vparameters.push_back("5.433,1.6506534,1.0841156,0.6295,0.2669,0.9708,0.2371,0.0205");  // 001, ternary metal-oxide prototype (ICSD #60825)
        vparameters.push_back("5.7187,1.5767395,1.2260479,0.6273,0.7275,0.0618,0.746,0.0217");  // 002, ternary metal-oxide prototype (ICSD #416147)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_oC28_63_c_acf_f"){
        vparameters.push_back("3.3255,3.3143888,3.1511051,0.8924,0.539,0.287,0.372,0.3652,0.9275");  // 001, ternary metal-oxide prototype (ICSD #79371)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oC28_63_f_c_acf"){
        vparameters.push_back("2.803,3.30396,3.3742419,0.1107,0.4523,0.6345,0.9267,0.73,0.1102");  // 001, ternary metal-oxide prototype (ICSD #161057)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oC24_64_d_f_f"){
        vparameters.push_back("8.618,0.74297981,0.74460432,0.1777,0.874,0.3746,0.8308,0.6668");  // 001, ternary metal-oxide prototype (ICSD #30964)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_oC28_64_a_f_ef"){
        vparameters.push_back("5.3565,2.4534117,1.0095398,0.5,0.3613,0.9912,0.1852,0.9939");  // 001, ternary metal-oxide prototype (ICSD #67836)
        vparameters.push_back("5.1374,2.4234048,1.0085841,0.48,0.13,0.49,0.34,0.51");  // 002, ternary metal-oxide prototype (ICSD #86754)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_oC28_64_a_ef_f"){
        vparameters.push_back("5.25,2.4281905,1.0131429,0.582,0.329,0.5,0.134,0.516");  // 001, ternary metal-oxide prototype (ICSD #91073)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC10_65_a_c_bf"){
        vparameters.push_back("5.4921,1.0131826,0.7048852");  // 001, ternary metal-oxide prototype (ICSD #28564)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC12_65_a_hi_c"){
        vparameters.push_back("6.5156,1.0347167,0.44390386,0.783,0.704");  // 001, ternary metal-oxide prototype (ICSD #72872)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC12_65_i_ai_c"){
        vparameters.push_back("3.9313,2.9412408,0.88840841,0.8342,0.6699");  // 001, ternary metal-oxide prototype (ICSD #99041)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_oC12_65_j_bj_a"){
        vparameters.push_back("4.4512,2.1387042,0.78237329,0.3444,0.1443");  // 001, ternary metal-oxide prototype (ICSD #154704) and Li2O3Pr (http://aflow.org/prototype-encyclopedia/A2B3C_oC12_65_h_bh_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oC16_65_b_ai_p"){
        vparameters.push_back("5.5,1.9963636,0.71272727,0.2521,0.2543,0.6264");  // 001, ternary metal-oxide prototype (ICSD #65881)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oC16_65_i_n_ad"){
        vparameters.push_back("4.838,1.6839603,1.2622985,0.25,0.662,0.254");  // 001, ternary metal-oxide prototype (ICSD #9389)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oC16_65_c_ai_p"){
        vparameters.push_back("5.47,2,0.71480804,0.2526,0.754,0.874");  // 001, ternary metal-oxide prototype (ICSD #83079)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C4_oC16_65_d_ai_p"){
        vparameters.push_back("5.54,1.9801444,0.71119134,0.248,0.248,0.627");  // 001, ternary metal-oxide prototype (ICSD #89232)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C3_oC20_65_b_jp_af"){
        vparameters.push_back("7.152,1.4045022,0.44010067,0.344,0.21,0.128");  // 001, ternary metal-oxide prototype (ICSD #35337)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5C2_oC20_65_ai_b2i_j"){
        vparameters.push_back("3.9375,4.9356698,0.88010159,0.8013,0.8989,0.2969,0.3995");  // 001, ternary metal-oxide prototype (ICSD #50089)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oC20_65_i_acm_j"){
        vparameters.push_back("5.291,1.8934039,0.99017199,0.8113,0.6737,0.748");  // 001, ternary metal-oxide prototype (ICSD #15927)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7C_oC22_65_af_cjp_b"){
        vparameters.push_back("7.904,1.4072621,0.48500759,0.135,0.32,0.12");  // 001, ternary metal-oxide prototype (ICSD #168909)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oC28_65_ij_ac_ijo"){
        vparameters.push_back("5.4687,2.3182841,0.99935999,0.63883,0.8225,0.13883,0.3225,0.25,0.746");  // 001, ternary metal-oxide prototype (ICSD #78338)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_67_m_g_a"){
        vparameters.push_back("6.5907,2.3255496,1.0527713,0.3111,0.08701,0.80131");  // 001, ternary metal-oxide prototype (ICSD #411905)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C2_oF28_69_h_ah_g"){
        vparameters.push_back("13.376,0.85690789,0.29388457,0.6288,0.8341,0.33");  // 001, ternary metal-oxide prototype (ICSD #68676)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_oF28_69_a_cg_g"){
        vparameters.push_back("12.8464,0.45590204,0.45093567,0.169,0.3591");  // 001, ternary metal-oxide prototype (ICSD #81577)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_oF48_70_h_b_a"){
        vparameters.push_back("10.239,0.61679852,0.97949018,0.2838,0.7006,0.2819");  // 001, ternary metal-oxide prototype (ICSD #88369)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oF56_70_d_a_h"){
        vparameters.push_back("8.3041,1.0022519,1.0058405,0.7613,0.26135,0.76093");  // 001, ternary metal-oxide prototype (ICSD #290589)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_oF56_70_a_g_h"){
        vparameters.push_back("12.783,0.50456857,0.84925291,0.794,0.3046,0.478,0.97");  // 001, ternary metal-oxide prototype (ICSD #151971)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oI10_71_e_f_a"){
        vparameters.push_back("9.354,0.39982895,0.31804576,0.735,0.357");  // 001, ternary metal-oxide prototype (ICSD #19007)
        vparameters.push_back("9.588,0.47820192,0.32530246,0.32,0.36");  // 002, ternary metal-oxide prototype (ICSD #25018)
        vparameters.push_back("9.8977,0.46022813,0.30390899,0.64242,0.6373");  // 003, ternary metal-oxide prototype (ICSD #51498)
        vparameters.push_back("9.3158,0.40290689,0.32007986,0.2997,0.8514");  // 004, ternary metal-oxide prototype (ICSD #61199)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_oI10_71_a_e_f"){
        vparameters.push_back("9.396,0.38984674,0.30470413,0.293,0.3565");  // 001, ternary metal-oxide prototype (ICSD #25001)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_oI12_71_b_e_af"){
        vparameters.push_back("10.967,0.39983587,0.32388073,0.8117,0.797");  // 001, ternary metal-oxide prototype (ICSD #10038)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI16_71_e_e_2f"){
        vparameters.push_back("5.757,0.79416363,0.48723293,0.876,0.3847,0.894,0.352");  // 001, ternary metal-oxide prototype (ICSD #15768)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C3_oI18_71_f_n_ae"){
        vparameters.push_back("13.175,0.41487666,0.23096774,0.3944,0.7025,0.1162,0.2557");  // 001, ternary metal-oxide prototype (ICSD #6157)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5C2_oI20_71_ae_def_e"){
        vparameters.push_back("18.662,0.21117244,0.18193656,0.7893,0.8948,0.60531,0.7831");  // 001, ternary metal-oxide prototype (ICSD #416904)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C7_oI28_71_aef_g_dhn"){
        vparameters.push_back("16.594,0.5605038,0.34175003,0.3482,0.6929,0.216,0.7888,0.1661,0.2219");  // 001, ternary metal-oxide prototype (ICSD #91309)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_oI24_72_c_bf_j"){
        vparameters.push_back("10.9,0.50091743,0.54311927,0.24,0.865,0.73");  // 001, ternary metal-oxide prototype (ICSD #16919)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_oI24_72_ce_b_j"){
        vparameters.push_back("9.945,0.60070387,0.57254902,0.0945,0.3268");  // 001, ternary metal-oxide prototype (ICSD #4204)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI32_72_j_f_2j"){
        vparameters.push_back("17.978,0.28863055,0.27689398,0.25,0.0893,0.2311,0.211,0.261,0.088,0.769");  // 001, ternary metal-oxide prototype (ICSD #25385)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oI16_74_f_e_c"){
        vparameters.push_back("9.46,0.57748414,0.49260042,0.7253,0.854");  // 001, ternary metal-oxide prototype (ICSD #2277)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oI20_74_e_ef_c"){
        vparameters.push_back("8.476,0.70870694,0.71531383,0.7467,0.1953,0.7187");  // 001, ternary metal-oxide prototype (ICSD #15933)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oI20_74_ac_i_e"){
        vparameters.push_back("6.1838,0.90243863,1.4640674,0.64160065,0.21467637,0.047724499");  // 001, ternary metal-oxide prototype (ICSD #36535)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_oI24_74_e_e_hi"){
        vparameters.push_back("7.173,0.99177471,0.87383243,0.375,0.875,0.554,0.812,0.305,0.965");  // 001, ternary metal-oxide prototype (ICSD #93789)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI28_79_c_2a_2c"){
        vparameters.push_back("8.484,0.68517209,0.1763,0.671,0.3317,0.1684,0.5,0.402,0.297,0.154,0.798,0.899,0.161");  // 001, ternary metal-oxide prototype (ICSD #12104)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI24_82_g_g_g"){
        vparameters.push_back("9.91,0.55095863,0.353,0.353,0.5,0.187,0.501,0.004,0.191,0.502,0.499");  // 001, ternary metal-oxide prototype (ICSD #25744)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C_tP14_85_a_cg_c"){
        vparameters.push_back("6.6078,0.64538878,0.7327,0.3395,0.7034,0.4623,0.7597");  // 001, ternary metal-oxide prototype (ICSD #27315)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C4_tI18_87_a_h_h"){
        vparameters.push_back("7.167,0.65762523,0.8038,0.5941,0.9185,0.2526");  // 001, ternary metal-oxide prototype (ICSD #67826)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI24_87_h_h_h"){
        vparameters.push_back("8.87,0.5264938,0.35,0.353,0.68,0,0.18,0.5");  // 001, ternary metal-oxide prototype (ICSD #15099)
        vparameters.push_back("9.56,0.60774059,0.35,0.353,0.18,0.5,0.68,0");  // 002, ternary metal-oxide prototype (ICSD #15100)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C4_tI26_87_a_2h_h"){
        vparameters.push_back("9.885,0.31633789,0.3015,0.6478,0.3422,0.0458,0.3323,0.8502");  // 001, ternary metal-oxide prototype (ICSD #1562)
        vparameters.push_back("9.963,0.29268293,0.0392,0.3349,0.6552,0.3015,0.8516,0.3307");  // 002, ternary metal-oxide prototype (ICSD #100596)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_tI32_87_h_di_e"){
        vparameters.push_back("8.935,0.67274762,0.81212,0.60628,0.81846,0.81385,0.90483,0.74225");  // 001, ternary metal-oxide prototype (ICSD #81)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_tI24_88_f_a_b"){
        vparameters.push_back("5.501,2.2035993,0.151,0.529,0.639");  // 001, ternary metal-oxide prototype (ICSD #16189)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tI24_88_a_f_b"){
        vparameters.push_back("5.47,2.226691,0.247,0.658,0.96");  // 001, ternary metal-oxide prototype (ICSD #39137)
        vparameters.push_back("5.64,2.251773,0.21,0.08,0.57");  // 002, ternary metal-oxide prototype (ICSD #52384)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tI24_88_b_f_a"){
        vparameters.push_back("5.611,2.2613616,0.2415,0.9914,0.7874");  // 001, ternary metal-oxide prototype (ICSD #169094)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tI28_88_f_d_a"){
        vparameters.push_back("5.8382,1.7417526,0.8478,0.9471,0.6755");  // 001, ternary metal-oxide prototype (ICSD #78843)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI28_88_d_a_f"){
        vparameters.push_back("6.4297,1.5943201,0.7665,0.3976,0.9218");  // 001, ternary metal-oxide prototype (ICSD #80327)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP16_92_a_a_b"){
        vparameters.push_back("5.1687,1.2126647,0.1759,0.8126,0.3369,0.2906,0.7723");  // 001, ternary metal-oxide prototype (ICSD #23815)
        vparameters.push_back("5.6208,1.3103651,0.1898,0.779,0.2293,0.8763,0.0455");  // 002, ternary metal-oxide prototype (ICSD #33768)
        vparameters.push_back("5.502,1.4187568,0.257,0.757,0.191,0.008,0.12");  // 003, ternary metal-oxide prototype (ICSD #151883)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tP5_99_ac_b_a"){
        vparameters.push_back("3.8072,1.2340303,0.2102,0.5674,0,0.6915");  // 001, ternary metal-oxide prototype (ICSD #152276)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tP5_99_b_a_ac"){
        vparameters.push_back("3.64,1.3200549,0.487,0.872,0.051,0.316");  // 001, ternary metal-oxide prototype (ICSD #188467)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B5C_tI24_107_bc_ad_a"){
        vparameters.push_back("7.707,0.7194758,0.55,0.94,0.75,0.207,0.51,0.728,0");  // 001, ternary metal-oxide prototype (ICSD #15102)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_tI28_108_c_d_a"){
        vparameters.push_back("8.623,0.68526035,0.8301,0.1778,0,0.402,0.29,0.35");  // 001, ternary metal-oxide prototype (ICSD #9622)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_113_c_a_e"){
        vparameters.push_back("5.5068,0.60928307,0.5627,0.2275,0.0917");  // 001, ternary metal-oxide prototype (ICSD #246244)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI24_119_g_i_i"){
        vparameters.push_back("8.515,0.44744568,0.842,0.69,0.553,0.659,0.04");  // 001, ternary metal-oxide prototype (ICSD #282)
        vparameters.push_back("9.893,0.55038916,0.352,0.678,0.0161,0.7058,0.5038");  // 002, ternary metal-oxide prototype (ICSD #24818)
        vparameters.push_back("9.52,0.48308824,0.8447,0.3385,0.4885,0.694,0.9893");  // 003, ternary metal-oxide prototype (ICSD #49752)
        vparameters.push_back("10.243,0.60216733,0.1417,0.6935,0.4971,0.2837,0.012");  // 004, ternary metal-oxide prototype (ICSD #49754)
        vparameters.push_back("8.514,0.44738078,0.8429,0.3592,0.4574,0.6941,0.9997");  // 005, ternary metal-oxide prototype (ICSD #49755)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_tI16_121_bd_i_a"){
        vparameters.push_back("5.93,1.3676223,0.8297,0.905");  // 001, ternary metal-oxide prototype (ICSD #4138)
        vparameters.push_back("4.9968,1.9394613,0.8073,0.9003");  // 002, ternary metal-oxide prototype (ICSD #417470)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_tI18_121_e_a_di"){
        vparameters.push_back("4.089,3.9104916,0.326,0.75,0.44");  // 001, ternary metal-oxide prototype (ICSD #25611)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C8_tI24_121_a_bd_2i"){
        vparameters.push_back("6.703,1.1385947,0.3645,0.6788,0.7921,0.0082");  // 001, ternary metal-oxide prototype (ICSD #9356) //DX20210428 - equivalent to K3CrO8 (http://aflow.org/prototype-encyclopedia/AB3C8_tI24_121_a_bd_2i.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC8_tI24_121_bd_a_2i"){
        vparameters.push_back("6.79,1.1605302,0.358,0.3,0.781,0.944");  // 001, ternary metal-oxide prototype (ICSD #30405)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_122_b_a_d"){
        vparameters.push_back("5.71,1.2767075,0.63");  // 001, ternary metal-oxide prototype (ICSD #4199)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI28_122_d_a_e"){
        vparameters.push_back("6.04,1.2880795,0.525,0.72,0.972,0.12");  // 001, ternary metal-oxide prototype (ICSD #16708)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_tP5_123_c_ae_b"){
        vparameters.push_back("3.982,1.0087896");  // 001, ternary metal-oxide prototype (ICSD #27965)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_tP5_123_a_be_d"){
        vparameters.push_back("3.993,1.010268");  // 001, ternary metal-oxide prototype (ICSD #109327)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tP5_123_c_b_ae"){
        vparameters.push_back("4.6322,0.50066923");  // 001, ternary metal-oxide prototype (ICSD #168904)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C6_tP11_123_a_eh_bci"){
        vparameters.push_back("4.169,1.9659391,0.2467,0.7558");  // 001, ternary metal-oxide prototype (ICSD #42006)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B9C2_tP16_123_bfg_cegi_h"){
        vparameters.push_back("4.166,2.9320691,0.8351,0.6598,0.6661,0.8272");  // 001, ternary metal-oxide prototype (ICSD #42007)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_tP10_127_c_bg_a"){
        vparameters.push_back("5.5552,0.70812932,0.7705");  // 001, ternary metal-oxide prototype (ICSD #23322)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_tP14_131_c_n_j"){
        vparameters.push_back("5.7786,0.9689025,0.242,0.268");  // 001, ternary metal-oxide prototype (ICSD #6159)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_tP12_136_ad_b_f"){
        vparameters.push_back("6.031,1.1865362,0.7054");  // 001, ternary metal-oxide prototype (ICSD #262579)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tI10_139_a_e_d"){
        vparameters.push_back("4.034,3.405057,0.344");  // 001, ternary metal-oxide prototype (ICSD #9099) and Co2S2Tl (http://aflow.org/prototype-encyclopedia/A2B2C_tI10_139_d_e_a.TlCo2S2.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tI10_139_a_e_e"){
        vparameters.push_back("3.42,3.8947368,0.675,0.853");  // 001, ternary metal-oxide prototype (ICSD #25511)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_tI14_139_a_e_cd"){
        vparameters.push_back("3.945,3.0851711,0.6487");  // 001, ternary metal-oxide prototype (ICSD #4203)
        vparameters.push_back("3.8329,3.4733492,0.6395");  // 002, ternary metal-oxide prototype (ICSD #261660)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_tI14_139_e_ce_a"){
        vparameters.push_back("4.296,3.0959032,0.645,0.845");  // 001, ternary metal-oxide prototype (ICSD #27113)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C6_tI22_139_ae_e_dg"){
        vparameters.push_back("3.9686,4.8670564,0.817,0.5826,0.5838");  // 001, ternary metal-oxide prototype (ICSD #249209)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI24_139_h_j_i"){
        vparameters.push_back("9.52,0.48497899,0.8442,0.6932,0.8405");  // 001, ternary metal-oxide prototype (ICSD #40153)
        vparameters.push_back("9.925,0.54992443,0.8522,0.2944,0.1772");  // 002, ternary metal-oxide prototype (ICSD #40154)
        vparameters.push_back("8.514,0.44738078,0.8431,0.3071,0.8624");  // 003, ternary metal-oxide prototype (ICSD #40156)
        vparameters.push_back("10.243,0.60216733,0.8582,0.713,0.8704");  // 004, ternary metal-oxide prototype (ICSD #40160)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B3C2_tI24_139_aeg_be_e"){
        vparameters.push_back("3.899,5.2321108,0.81,0.685,0.905,0.405");  // 001, ternary metal-oxide prototype (ICSD #20294) //DX20210428 - equivalent to O7Sr3Ti2 (http://aflow.org/prototype-encyclopedia/A7B3C2_tI24_139_aeg_be_e.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_tI24_139_h_i_j"){
        vparameters.push_back("9.541,0.60916047,0.1354,0.2697,0.1923");  // 001, ternary metal-oxide prototype (ICSD #40159)
        vparameters.push_back("9.937,0.5856194,0.8518,0.622,0.818");  // 002, ternary metal-oxide prototype (ICSD #409552)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C8_tI30_139_2e_ae_cdg"){
        vparameters.push_back("3.9168,6.4649714,0.5658,0.7012,0.8746,0.8733");  // 001, ternary metal-oxide prototype (ICSD #51097)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tI20_140_ah_b_c"){
        vparameters.push_back("5.8504,1.417886,0.2222");  // 001, ternary metal-oxide prototype (ICSD #1522)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_tI20_140_c_ah_b"){
        vparameters.push_back("8.322,1,0.303");  // 001, ternary metal-oxide prototype (ICSD #164621)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI16_141_b_a_e"){
        vparameters.push_back("4.66,2.2575107,0.1529");  // 001, ternary metal-oxide prototype (ICSD #22040)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tI20_141_a_d_e"){
        vparameters.push_back("5.722,1.7588256,0.135");  // 001, ternary metal-oxide prototype (ICSD #9456)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tI24_141_a_h_b"){
        vparameters.push_back("5.9607,2.4235409,0.013,0.3114");  // 001, ternary metal-oxide prototype (ICSD #72817)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3_tI28_141_d_c_ae"){
        vparameters.push_back("5.8978,1.8166096,0.1131");  // 001, ternary metal-oxide prototype (ICSD #51502)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_tI28_141_d_h_a"){
        vparameters.push_back("5.722,1.6141209,0.523,0.739");  // 001, ternary metal-oxide prototype (ICSD #15305)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_tI28_141_a_d_h"){
        vparameters.push_back("5.81,1.6987952,0.55,0.725");  // 001, ternary metal-oxide prototype (ICSD #24258)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tI28_141_d_a_h"){
        vparameters.push_back("5.76,1.4756944,0.989,0.733");  // 001, ternary metal-oxide prototype (ICSD #37023)
        vparameters.push_back("5.87,1.3475298,0.496,0.28");  // 002, ternary metal-oxide prototype (ICSD #183968)
        vparameters.push_back("6.0136,1.2943495,0.9656,0.2555");  // 003, ternary metal-oxide prototype (ICSD #246080)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tI32_141_c_d_h"){
        vparameters.push_back("5.6504,1.6356364,0.983,0.752");  // 001, ternary metal-oxide prototype (ICSD #40486)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI32_141_d_h_c"){
        vparameters.push_back("4.043,2.1340589,0.5,0.2362");  // 001, ternary metal-oxide prototype (ICSD #164158)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR5_146_a_a_b"){
        vparameters.push_back("5.4886776,2.4295959,0,0.231,0.943,0.398,0.542");  // 001, ternary metal-oxide prototype (ICSD #162264)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B9C4_hR16_146_3a_3b_4a"){
        vparameters.push_back("6.0420286,4.1143665,0,0.16914,0.58535,0.44044,0.87605,0.74433,0.31543,0.0343,0.569,0.18,0.6202,0.7103,0.1201,0.4162,0.6314,0.079");  // 001, ternary metal-oxide prototype (ICSD #33239)
        vparameters.push_back("6.1100281,4.1222577,0,0.836,0.421,0.573,0.142,0.265,0.699,0.9,0.94,0.41,0.452,0.882,0.252,0.555,0.955,0.415");  // 002, ternary metal-oxide prototype (ICSD #87118)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C9_hR16_146_3a_4a_3b"){
        vparameters.push_back("6.0975041,4.1221864,0,0.8383,0.4251,0.5755,0.1351,0.2629,0.7042,0.858,0.366,0.011,0.837,0.323,0.411,0.978,0.382,0.589");  // 001, ternary metal-oxide prototype (ICSD #33807)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_hP11_147_a_dg_d"){
        vparameters.push_back("5.8579,1.03561,0.57,0.2891,0.135,0.349,0.206");  // 001, ternary metal-oxide prototype (ICSD #51014)
        vparameters.push_back("5.6669,1.082726,0.602,0.2799,0.12,0.325,0.186");  // 002, ternary metal-oxide prototype (ICSD #51016)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B6C_hR15_148_cf_f_a"){
        vparameters.push_back("5.461035,2.7976535,0.6574,0.8996,0.5006,0.2356,0.2341,0.6161,0.9071");  // 001, ternary metal-oxide prototype (ICSD #1180) //DX20210428 - equivalent to Li7O6Ta (partially occupied Li, http://aflow.org/prototype-encyclopedia/A8B6C_hR15_148_cf_f_a.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP12_152_b_c_a"){
        vparameters.push_back("5.886,1.1440707,0.582,0.653,0.118,0.42,0.87335");  // 001, ternary metal-oxide prototype (ICSD #25559)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP12_152_a_c_b"){
        vparameters.push_back("5.927,1.1316011,0.35,0.5,0.19,0.78,0.98336667");  // 001, ternary metal-oxide prototype (ICSD #25812)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP12_154_a_c_b"){
        vparameters.push_back("3.833,4.9543439,0.375,0.371,0.631,0.628,0.10296667");  // 001, ternary metal-oxide prototype (ICSD #69739)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR5_160_a_b_a"){
        vparameters.push_back("5.650857,1.2295576,0.013,0.5,0.524,0.031");  // 001, ternary metal-oxide prototype (ICSD #6102)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C2_hR7_160_2a_b_2a"){
        vparameters.push_back("6.000047,2.3883328,0.7445,0.238,0,0.482,0.2303,0.7238");  // 001, ternary metal-oxide prototype (ICSD #2216)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR10_161_a_b_a"){
        vparameters.push_back("5.442009,2.5242564,0,0.2226,0.1536,0.2924,0.8333");  // 001, ternary metal-oxide prototype (ICSD #187684)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_hR12_161_2a_b_a"){
        vparameters.push_back("4.9711028,2.9747927,0.688,0.831,0,0.3949,0.1029,0.7369");  // 001, ternary metal-oxide prototype (ICSD #200999)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B6C2_hP13_162_ef_k_d"){
        vparameters.push_back("5.9324,1.0805913,0.2413,0.62216667,0.6889");  // 001, ternary metal-oxide prototype (ICSD #40058)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C2_hP13_162_de_ak_c"){
        vparameters.push_back("5.628,1.3041933,0.726,0.64896667,0.662");  // 001, ternary metal-oxide prototype (ICSD #6321)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_hP6_164_b_e_d"){
        vparameters.push_back("5.25,1.0365714,0.19");  // 001, ternary metal-oxide prototype (ICSD #413342)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C_hP11_164_d_di_a"){
        vparameters.push_back("5.846,1.0253336,0.70817,0.42828,0.83386,0.19598");  // 001, ternary metal-oxide prototype (ICSD #59999)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B4C2_hP13_164_ai_2d_c"){
        vparameters.push_back("5.9388,1.3019802,0.7677,0.0837,0.5826,0.1568,0.2991");  // 001, ternary metal-oxide prototype (ICSD #72809)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8C3_hP14_164_ad_di_bd"){
        vparameters.push_back("5.651,1.3080871,0.312,0.957,0.739,0.163,0.658");  // 001, ternary metal-oxide prototype (ICSD #100782)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_hP14_164_di_abd_d"){
        vparameters.push_back("6.278,1.2900605,0.525,0.8283,0.3145,0.18096667,0.228");  // 001, ternary metal-oxide prototype (ICSD #8212)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR4_166_b_a_c"){
        vparameters.push_back("2.8627267,5.9309897,0.12");  // 001, ternary metal-oxide prototype (ICSD #25593)
        vparameters.push_back("2.9200073,5.8630223,0.245");  // 002, ternary metal-oxide prototype (ICSD #166514)
        vparameters.push_back("3.7008243,3.9759849,0.415");  // 003, ternary metal-carbo-nitride prototype (ICSD #25763)
        vparameters.push_back("3.8144175,3.6265766,0.37");  // 004, ternary metal-carbo-nitride prototype (ICSD #31100)
        vparameters.push_back("3.2734244,4.316063,0.4117");  // 005, ternary metal-carbo-nitride prototype (ICSD #75039)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hR4_166_c_a_b"){
        vparameters.push_back("3.2500058,5.9138415,0.782");  // 001, ternary metal-oxide prototype (ICSD #31959)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR5_166_b_a_d"){
        vparameters.push_back("5.5757966,1.2408123");  // 001, ternary metal-oxide prototype (ICSD #20288)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_hR5_166_c_a_c"){
        vparameters.push_back("2.9199816,8.2561836,0.21382,0.62405");  // 001, ternary metal-oxide prototype (ICSD #160574)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_hR7_166_c_2c_a"){
        vparameters.push_back("3.4550125,7.251532,0.785,0.7075,0.8708");  // 001, ternary metal-oxide prototype (ICSD #4192)
        vparameters.push_back("3.5135799,7.0520644,0.21414,0.292,0.1279");  // 002, ternary metal-oxide prototype (ICSD #67701)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C2_hR7_166_ab_d_c"){
        vparameters.push_back("6.0010193,2.3874335,0.24013");  // 001, ternary metal-oxide prototype (ICSD #15511) //DX20210428 - equivalent to Ni3Pb2S2 (shandite, http://aflow.org/prototype-encyclopedia/A3B2C2_hR7_166_d_ab_c.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_hR7_166_c_a_2c"){
        vparameters.push_back("3.3383925,7.8070463,0.21636,0.12799,0.29328");  // 001, ternary metal-oxide prototype (ICSD #157323)
        vparameters.push_back("4.309,9.44047342771,0.4247,0.1332,0.294");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hR8_166_c_ad_c"){
        vparameters.push_back("5.7093278,2.4542326,0.12499595,0.37721836");  // 001, ternary metal-oxide prototype (ICSD #29010)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_hR12_166_c_h_ae"){
        vparameters.push_back("6.0730183,3.0918812,0.1686,0.1224,0.6049");  // 001, ternary metal-oxide prototype (ICSD #248051)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C8_hR13_166_a_bd_ch"){
        vparameters.push_back("5.6699893,3.3068807,0.78,0.28,0.78");  // 001, ternary metal-oxide prototype (ICSD #20898)
        vparameters.push_back("5.6199759,2.916373,0.222,0.722,0.222");  // 002, ternary metal-oxide prototype (ICSD #40470)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C8_hR13_166_ac_c_ch"){
        vparameters.push_back("5.7570174,3.7151323,0.2053,0.408,0.327,0.738,0.218");  // 001, ternary metal-oxide prototype (ICSD #9457)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8C3_hR13_166_c_ch_e"){
        vparameters.push_back("6.1959757,2.4861206,0.145,0.397,0.886,0.44");  // 001, ternary metal-oxide prototype (ICSD #65412)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B8C_hR13_166_bd_ch_a"){
        vparameters.push_back("5.6599993,3.1254404,0.221,0.721,0.221");  // 001, ternary metal-oxide prototype (ICSD #40469)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B3C2_hR13_166_ch_ac_c"){
        vparameters.push_back("5.9818732,3.4233652,0.335,0.208,0.412,0.282,0.75");  // 001, ternary metal-oxide prototype (ICSD #27651)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_hR14_166_c_ch_ad"){
        vparameters.push_back("5.7530302,2.5072375,0.3725,0.7485,0.7553,0.245");  // 001, ternary metal-oxide prototype (ICSD #151457)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hR15_166_bc_dh_ac"){
        vparameters.push_back("5.7500172,3.7565182,0.2819,0.1181,0.786,0.254");  // 001, ternary metal-oxide prototype (ICSD #10253)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4_hR16_167_e_a_be"){
        vparameters.push_back("7.860865,2.1035757,0.6667,0");  // 001, ternary metal-oxide prototype (ICSD #30662)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_hP14_182_b_f_cg"){
        vparameters.push_back("5.458,1.6518871,0.055,0.343");  // 001, ternary metal-oxide prototype (ICSD #2769)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_186_b_a_c"){
        vparameters.push_back("5.631,0.85384479,0.0275,0.25,0.15003333,0.725");  // 001, ternary metal-oxide prototype (ICSD #15761)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hP10_186_a_a2b_b"){
        vparameters.push_back("3.909,3.1632131,0,0.2,0.07,0.35,0.755");  // 001, ternary metal-oxide prototype (ICSD #200088)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4_hP14_186_b_ab_bc"){
        vparameters.push_back("6.203,1.2248912,0,0.677,0.074,0.459,0.81193333,0.25");  // 001, ternary metal-oxide prototype (ICSD #39504)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_187_bc_gh_i"){
        vparameters.push_back("3.52,3.8636364,0.16,0.34,0.16");  // 001, ternary metal-oxide prototype (ICSD #49652)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_187_ad_i_jk"){
        vparameters.push_back("5.645,0.84180691,0.25594,0.5185,0.8149");  // 001, ternary metal-oxide prototype (ICSD #88670)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C_hP16_190_af_cg_d"){
        vparameters.push_back("4.5961,1.9429516,0.42925,0.33371667");  // 001, ternary metal-oxide prototype (ICSD #26179)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_hP6_191_d_f_a"){
        vparameters.push_back("5.321,0.67242999");  // 001, ternary metal-oxide prototype (ICSD #15995)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_hP8_194_f_c_a"){
        vparameters.push_back("3.251,3.9329437,0.5837");  // 001, ternary metal-oxide prototype (ICSD #1270)
        vparameters.push_back("3.2416,2.8226802,0.6624");  // 002, ternary metal-nitride prototype (ICSD #80373)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP8_194_c_f_a"){
        vparameters.push_back("3.5206,3.2431972,0.0893");  // 001, ternary metal-oxide prototype (ICSD #35580)
        vparameters.push_back("3.223,3.5411108,0.0893");  // 002, ternary metal-oxide prototype (ICSD #60847)
        vparameters.push_back("3.1533,2.9404116,0.383");  // 003, ternary metal-carbo-nitride prototype (ICSD #249388)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hP10_194_c_bf_a"){
        vparameters.push_back("3.68,2.8586957,0.077");  // 001, ternary metal-oxide prototype (ICSD #27100)
        vparameters.push_back("3.61,3.1551247,0.085");  // 002, ternary metal-oxide prototype (ICSD #73361)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_194_c_a_bf"){
        vparameters.push_back("3.31,3.6371601,0.0861");  // 001, ternary metal-oxide prototype (ICSD #30339)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_194_a_c_bf"){
        vparameters.push_back("3.3985,3.3765485,0.0871");  // 001, ternary metal-oxide prototype (ICSD #67671)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hP16_194_ac_f_ef"){
        vparameters.push_back("2.875,8.5213913,0.0858,0.625,0.1642");  // 001, ternary metal-oxide prototype (ICSD #261608)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C2_cI28_199_a_b_a"){
        vparameters.push_back("8.419,0.5232,0.7562,0.2887");  // 001, ternary metal-oxide prototype (ICSD #1412)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C7_cF52_216_e_bc_ag"){
        vparameters.push_back("8.09,0.8658,0.007");  // 001, ternary metal-oxide prototype (ICSD #100355)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C_cP16_223_e_c_a"){
        vparameters.push_back("5.826");  // ternary metal-oxide prototype (ICSD #16537)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_cF40_225_c_ab_e"){
        vparameters.push_back("8.7676,0.7395");  // 001, ternary metal-oxide prototype (ICSD #68714)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B6C_cF40_225_bc_e_a"){
        vparameters.push_back("8.309,0.75");  // 001, ternary metal-oxide prototype (ICSD #4148)
        vparameters.push_back("8.62,0.8");  // 002, ternary metal-oxide prototype (ICSD #76437)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC8_cF60_225_d_a_ce"){
        vparameters.push_back("8.381,0.77");  // 001, ternary metal-oxide prototype (ICSD #24710)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B8C_cF60_225_d_f_a"){
        vparameters.push_back("9.314,0.855");  // 001, ternary metal-oxide prototype (ICSD #280596)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C_cF64_225_f_d_ab"){
        vparameters.push_back("9.596,0.861");  // 001, ternary metal-oxide prototype (ICSD #2275)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_cF64_227_c_e_d"){
        vparameters.push_back("8.3756,0.2552");  // 001, ternary metal-oxide prototype (ICSD #48128)
      }
      // -------------------------------------------------------------------------
      // metal-carbide prototypes (from DX) //DX20210120
      // ---------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_aP12_2_4i_2i"){
        vparameters.push_back("8.5362,1.3554275,0.45471053,89.5283,98.7891,92.1312,0.064242394,0.12407192,0.3313133,0.18260225,0.1319026,0.17964091,0.6615567,0.089094821,0.19992363,0.58493164,0.16342196,0.32569874,0.38005957,0.13186639,0.74732262,0.86498513,0.10470324,0.76891375");  // 001, binary metal-carbide prototype (ICSD #66663)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4_mP10_11_e_4e"){
        vparameters.push_back("5.423,0.40494191,0.97602803,103.01,0.743,0.102,0.156,0.124,0.489,0.816,0.074,0.671,0.678,0.522");  // 001, binary metal-carbide prototype (ICSD #187143)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B6_mC22_12_agh_ij"){
        vparameters.push_back("5.4605,1.7320575,1,109.47,0.667,0.833,0.75,0.25,0.75,0.667,0.25");  // 001, binary metal-carbide prototype (ICSD #20695)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_mP10_14_2e_a"){
        vparameters.push_back("6.2877,0.58748032,0.88544301,140.6903,0.259,0.295,0.009,0.398,0.897,0.801");  // 001, binary metal-carbide prototype (ICSD #181493)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC12_15_f_e"){
        vparameters.push_back("6.5995,0.63488143,1.2562618,122.5269,0.5679,0.5242,0.604,0.0562");  // 001, binary metal-carbide prototype (ICSD #54184)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5_oP16_55_ah_d2g"){
        vparameters.push_back("4.9014,2.456543,0.63651202,0.612,0.348,0.67,0.107,0.848,0.28");  // 001, binary metal-carbide prototype (ICSD #181485)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oP10_58_ag_g"){
        vparameters.push_back("5.2786,1.2144887,0.7063047,0.207,0.1188,0.3901,0.7903");  // 001, binary metal-carbide prototype (ICSD #71941)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP16_62_2c_2c"){
        vparameters.push_back("9.0484,0.53451439,1.0976747,0.7102,0.0848,0.8249,0.1501,0.1591,0.0713,0.5084,0.8298");  // 001, binary metal-carbide prototype (ICSD #51529)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP16_62_3c_c"){
        vparameters.push_back("13.4651,0.19123512,0.3183712,0.01257,0.6471,0.38157,0.6804,0.29295,0.022,0.15034,0.8502");  // 001, binary metal-carbide prototype (ICSD #181492)
        vparameters.push_back("10.819,0.613735095665,0.89305850818,0.84315,0.15287,0.6188,-0.04894,0.42157,0.76482,0.16484,0.02957");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_oC20_68_i_a"){
        vparameters.push_back("3.9209,2.3473437,1.2146701,0.83398,0.37663,0.91596");  // 001, binary metal-carbide prototype (ICSD #181496)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_tP3_123_h_a"){
        vparameters.push_back("2.691,1.427722,0.7");  // 001, binary metal-carbide prototype (ICSD #184664)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_tI32_142_f_e"){
        vparameters.push_back("6.743,1.8795788,0.19,0.9371");  // 001, binary metal-carbide prototype (ICSD #28066)
        vparameters.push_back("7.58,1.9379947,0.785,0.063");  // 002, binary metal-carbide prototype (ICSD #36142)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_hP15_152_2ac_b"){
        vparameters.push_back("4.8258,1.2064114,0.077,0.61,0.318,0.489,0.221,0.55603333");  // 001, binary metal-carbide prototype (ICSD #181495)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_hR7_160_4a_3a"){
        vparameters.push_back("3.3301727,7.4772914,0.705,0.129,0.869,0.296,0,0.781,0.217");  // 001, binary metal-carbide prototype (ICSD #14397)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B8_hR13_160_2ab_2a2b"){
        vparameters.push_back("6.121,2.4494897,0,0.5,0.752,0.241,0.0001,0.5,0.7434,0.2433,0.2496,0.7497");  // 001, binary metal-carbide prototype (ICSD #618924)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B8_hR13_166_abd_ch"){
        vparameters.push_back("6.1152887,2.4363266,0.75383333,0.74656666,0.24683332");  // 001, binary metal-carbide prototype (ICSD #20822)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_187_a_i"){
        vparameters.push_back("3.0029,1.574811,0.25");  // 001, binary metal-carbide prototype (ICSD #43669)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP3_191_a_c"){
        vparameters.push_back("2.993,0.77246909");  // 001, binary metal-carbide prototype (ICSD #162103)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B_hP9_191_cl_b"){
        vparameters.push_back("4.945,1.1971689,0.16673");  // 001, binary metal-carbide prototype (ICSD #74641)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B_hP14_194_i_c"){
        vparameters.push_back("4.3,2.2651163,0.33336667");  // 001, binary metal-carbide prototype (ICSD #169041)
        vparameters.push_back("4.32,2.1173611,0.33336667");  // 002, binary metal-carbide prototype (ICSD #601565)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_cI28_220_a_c"){
        vparameters.push_back("7.2081,0.0496");  // 001, binary metal-carbide prototype (ICSD #42760)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cP12_224_e_b"){
        vparameters.push_back("5.74,0.65");  // 001, binary metal-carbide prototype (ICSD #31092)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF64_227_cd_e"){
        vparameters.push_back("8.634,0.248");  // 001, binary metal-carbide prototype (ICSD #618950)
      }
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      // ternaries
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_mC16_9_2a_a_a"){
        vparameters.push_back("5.327,0.99981228,1.5402478,117.8697,0.974,0.984,0.307,0.967,0.821,0.428,0.228,0.997,0.79,0.062,0.75,0.062");  // 001, ternary metal-nitride prototype (ICSD #55570)
        vparameters.push_back("5.3488,1.0005235,1.5376346,117.9499,0.061,0.676,0.578,0.056,0.562,0.729,0.8639,0.6796,0.3065,0.5,0.3223,0");  // 002, ternary metal-nitride prototype (ICSD #67375)
        vparameters.push_back("7.401,0.72395622,0.72382111,78.27,0.773,0.82,0.171,0.427,0.425,0.02,0.701,0.061,0.941,0.5,0.573,0");  // 003, ternary metal-nitride prototype (ICSD #658525)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B4C3_mC24_12_b2i_2i_ai"){
        vparameters.push_back("11.087,0.34012808,0.88743574,67.42,0.667,0.392,0.608,0.311,0.7442,0.8396,0.9006,0.3971,0.541,0.1638");  // 001, ternary metal-nitride prototype (ICSD #71440)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_mC32_12_2j_i_ghi"){
        vparameters.push_back("5.5319,2.171207,1.2202498,127.5955,0.81259,0.18844,0.99927,0.73337,0.51334,0.75735,0.3345,0.12551,0.914,0.3346,0.87524,0.4142");  // 001, ternary metal-nitride prototype (ICSD #236391)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_mC32_12_2j_i_bc2i"){
        vparameters.push_back("12.814,0.35008584,0.99996878,150.0273,0.73824,0.73824,0.09867,0.9095,0.59923,0.4083,0.8137,0.34,0.685,0.6852,0.3413,0.815");  // 001, ternary metal-nitride prototype (ICSD #420074)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C4_oP11_25_c2g_g_adh"){
        vparameters.push_back("3.5288,3.2390614,1.0376332,0,0.508,0.7921,0.2285,0.999,0.331,0.802,0.38511,0.301,0.20241,0.4972");  // 001, ternary metal-nitride prototype (ICSD #73156)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC8_38_d_b_a"){
        vparameters.push_back("3.648,1.2878289,1.8083882,0.632,0,0.67,0.86");  // 001, ternary metal-nitride prototype (ICSD #63064)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP16_57_2d_d_c"){
        vparameters.push_back("8.0662,0.83435819,0.80325308,0.03,0.752,0.557,0.877,0.46,0.346,0.5351");  // 001, ternary metal-nitride prototype (ICSD #189824)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP16_59_2e_ab_e"){
        vparameters.push_back("3.2434,2.8211136,1.7311463,0.072,0.4057,0.8525,0.082,0.3279,0.417,0,0.749");  // 001, ternary metal-nitride prototype (ICSD #80372)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oC16_63_g_c_a"){
        vparameters.push_back("4.365,2.125315,1.1917526,0.31705,0.658,0.087");  // 001, ternary metal-nitride prototype (ICSD #80312)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oC20_63_ac_f_c"){
        vparameters.push_back("2.854,3.2529783,2.4523476,0.745,0.094,0.642,0.931");  // 001, ternary metal-nitride prototype (ICSD #20297)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_oC16_65_p_a_bi"){
        vparameters.push_back("4.4794,2.7586284,0.37138903,0.80904,0.1614,0.62861");  // 001, ternary metal-nitride prototype (ICSD #420076)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C5_oC28_65_aegj_i_bq"){
        vparameters.push_back("7.7927,1.7478024,0.41143634,0.82651,0.283812,0.38165,0.28909,0.628546");  // 001, ternary metal-nitride prototype (ICSD #72382)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC3_oI16_71_n_b_af"){
        vparameters.push_back("12.101,0.36213536,0.27847285,0.18893,0.37493,0.6656");  // 001, ternary metal-nitride prototype (ICSD #72863)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_oI28_72_2j_j_a"){
        vparameters.push_back("9.4262,0.7963018,0.53107297,0.8689,0.3422,0.9219,0.1772,0.6459,0.1562");  // 001, ternary metal-nitride prototype (ICSD #42968)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP4_123_g_a_d"){
        vparameters.push_back("4.4612,1.102394,0.378");  // 001, ternary metal-nitride prototype (ICSD #391118)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tP4_123_a_g_d"){
        vparameters.push_back("4.2267,1.2508576,0.3843");  // 001, ternary metal-nitride prototype (ICSD #410874)
        vparameters.push_back("3.73,1.3941019,0.3846");  // 002, ternary metal-nitride prototype (ICSD #411254)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_tP5_123_b_e_ac"){
        vparameters.push_back("3.9406,0.97624727");  // 001, ternary metal-nitride prototype (ICSD #44353)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP8_129_2c_a_c"){
        vparameters.push_back("3.3344,2.1869002,0.175,0.348,0.653");  // 001, ternary metal-nitride prototype (ICSD #62598)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tP8_131_a_j_f"){
        vparameters.push_back("5.2467,1.6254026,0.6159");  // 001, ternary metal-nitride prototype (ICSD #410873)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP8_131_j_a_f"){
        vparameters.push_back("4.9278,1.5277609,0.6219");  // 001, ternary metal-nitride prototype (ICSD #412038)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5C2_hR11_160_4a_5a_2a"){
        vparameters.push_back("3.3216992,12.327423,0.9304,0.8248,0.648,0.5498,0,0.8717,0.734,0.6007,0.4697,0.3702,0.1016");  // 001, ternary metal-nitride prototype (ICSD #173676)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C3_hR13_160_4a_6a_3a"){
        vparameters.push_back("3.3180704,14.760707,0.0731,0.1545,0.4474,0.5362,0,0.112,0.2194,0.3871,0.4976,0.6091,0.3047,0.6925,0.9173");  // 001, ternary metal-nitride prototype (ICSD #173677)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_hP5_164_c_d_a"){
        vparameters.push_back("4.464,1.1798835,0.384,0.282");  // 001, ternary metal-nitride prototype (ICSD #50172)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_hP5_164_c_a_d"){
        vparameters.push_back("5.3581,0.98706631,0.3809,0.2775");  // 001, ternary metal-nitride prototype (ICSD #94394)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5C2_hR11_166_2c_a2c_c"){
        vparameters.push_back("3.2740153,12.28464,0.0891,0.8116,0.7326,0.1359,0.3668");  // 001, ternary metal-nitride prototype (ICSD #161586)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B6C3_hR13_166_2c_3c_ac"){
        vparameters.push_back("3.3020195,14.76676,0.76756667,0.1476,0.69396667,0.0823,0.80746667,0.38746667");  // 001, ternary metal-nitride prototype (ICSD #161585)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_hP14_186_3b_3b_a"){
        vparameters.push_back("3.355,5.0002981,0,0.6297,0.8706,0.2582,0.7504,0.0848,0.4103");  // 001, ternary metal-nitride prototype (ICSD #43477)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_hP4_187_a_g_d"){
        vparameters.push_back("3.7882,1.4064727,0.38");  // 001, ternary metal-nitride prototype (ICSD #410868)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_hP12_194_b_f_af"){
        vparameters.push_back("3.072,6.0970052,0.0701,0.629");  // 001, ternary metal-nitride prototype (ICSD #93503)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_hP12_194_f_b_af"){
        vparameters.push_back("3.12,5.9403846,0.0671,0.1234");  // 001, ternary metal-nitride prototype (ICSD #160572)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C_hP14_194_cf_df_a"){
        vparameters.push_back("3.355,5.0002981,0.6293,0.0876");  // 001, ternary metal-nitride prototype (ICSD #62308)
      }
      // -------------------------------------------------------------------------
      // carbo-nitride prototypes (from DX) //DX20210219
      // ---------------------------------------------------------------------------
      // binaries 
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B4_oP15_16_acgqtu_u"){
        vparameters.push_back("3.5013,1.029332,2.0983635,0.2793,0.2557,0.2324,0.2354,0.1394,0.7553,0.2671,0.381");  // 001, binary carbo-nitride prototype (ICSD #184897)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_oC12_36_a_2a"){
        vparameters.push_back("2.481,2.5961306,1.8476421,0.11253,0.22046,0.07558,0.5479,0.26266,0.67879");  // 001, binary carbo-nitride prototype (ICSD #247680)
        vparameters.push_back("3.145,3.55612082671,1.84702384738,0.0934,0.0,0.294,0.803,0.426,0.121");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B4_tP15_111_abcmn_n"){
        vparameters.push_back("3.5421,1.9760312,0.7572,0.2552,0.6281,0.2551,0.126");  // 001, binary carbo-nitride prototype (ICSD #184896)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tP6_113_a_e"){
        vparameters.push_back("3.313,1.4331422,0.34558,0.8392");  // 001, binary carbo-nitride prototype (ICSD #247676)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI6_119_a_f"){
        vparameters.push_back("2.469,2.4264885,0.86191");  // 001, binary carbo-nitride prototype (ICSD #247678)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI24_122_d_e"){
        vparameters.push_back("6.476,0.5302656,0.55701,0.43086,0.18919,0.46428");  // 001, binary carbo-nitride prototype (ICSD #247677)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_hR7_160_b_ab"){
        vparameters.push_back("4.7441335,1.9377301,0.665,0.491,0.018,0.17,0.66");  // 001, binary carbo-nitride prototype (ICSD #41952)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_164_c_2d"){
        vparameters.push_back("2.363,4.5222175,0.07383,0.46686,0.87827");  // 001, binary carbo-nitride prototype (ICSD #247679)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_hP14_176_h_ch"){
        vparameters.push_back("6.41,0.375039,0.178,0.772,0.331,0.033");  // 001, binary carbo-nitride prototype (ICSD #41950) //DX20210428 - equivalent to beta-Si3N4 (http://aflow.org/prototype-encyclopedia/A4B3_hP14_176_ch_h.html, part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_hP14_187_jk_adjk"){
        vparameters.push_back("4.742,1.417229,0.1759,0.4974,0.5099,0.8307");  // 001, binary carbo-nitride prototype (ICSD #83265)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_cP7_215_c_e"){
        vparameters.push_back("3.43,0.254");  // 001, binary carbo-nitride prototype (ICSD #41951)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_cI28_220_b_c"){
        vparameters.push_back("5.3973,0.7812");  // 001, binary carbo-nitride prototype (ICSD #83263)
      }
      // ---------------------------------------------------------------------------
      // ---------------------------------------------------------------------------
      // ternaries 
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_aP15_2_ai_3i_3i"){
        vparameters.push_back("6.6255,1.4557392,0.80572032,81.315,113.178,81.235,0.738,0.358,0.276,0.007,0.926,0.213,0.712,0.414,0.452,0.246,0.702,0.898,0.5673,0.68249,0.3839,0.9465,0.68114,0.0841,0.6833,0.04635,0.3214");  // 001, ternary metal-carbo-nitride prototype (ICSD #417297)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3C_aP16_2_4i_3i_i"){
        vparameters.push_back("3.705,2.4825911,2.0893387,105.76,100.16,88.31,0.714,0.7073,0.2955,0.7646,0.5506,0.2379,0.6416,0.7724,0.4729,0.7492,0.7996,0.1807,0.8074,0.423,0.1945,0.5842,0.8236,0.6184,0.7789,0.8762,0.0863,0.8207,0.1596,0.1566");  // 001, ternary metal-carbo-nitride prototype (ICSD #31926)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mC6_8_a_a_a"){
        vparameters.push_back("5.5307,0.94188801,0.67678594,85.58,0.574,0.635,0,0,0.417,0.456");  // 001, ternary metal-carbo-nitride prototype (ICSD #27350)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mC12_9_a_a_a"){
        vparameters.push_back("9.103,0.50598704,0.89220037,129.2616,0.286,0.836,0.429,0,0.832,0,0.741,0.168,0.495");  // 001, ternary metal-carbo-nitride prototype (ICSD #173942)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP16_33_a_2a_a"){
        vparameters.push_back("5.5531,2.1127478,0.69640381,0.3792,0.8977,0.462,0.1732,0.9407,0.452,0.56987,0.8555,0.47525,0.89288,0.86409,0");  // 001, ternary metal-carbo-nitride prototype (ICSD #16600)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C3_oI32_46_b_2bc_bc"){
        vparameters.push_back("7.98,1.2769424,0.77819549,0.1519,0,0.943,0.387,0.87,0.604,0.007,0.241,0.094,0.655,0.204,0.964,0.686,0.271");  // 001, ternary metal-carbo-nitride prototype (ICSD #43823)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP6_59_a_a_b"){
        vparameters.push_back("3.63,1.5013774,1.3360882,0.864,0.626,0.75");  // 001, ternary metal-carbo-nitride prototype (ICSD #77172)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oC16_63_c_a_f"){
        vparameters.push_back("2.9916,2.0551544,3.1441369,0.6153,0.3802,0.619");  // 001, ternary metal-carbo-nitride prototype (ICSD #161460)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI32_122_d_e_d"){
        vparameters.push_back("8.8047,0.61704544,0.0675,0.54692,0.8084,0.428,0.5795");  // 001, ternary metal-carbo-nitride prototype (ICSD #280523)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C6_hR11_155_d_c_f"){
        vparameters.push_back("6.2732404,2.3402753,0.8284,0.7025,0.9907,0.5925,0.1807");  // 001, ternary metal-carbo-nitride prototype (ICSD #240311)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_hR4_160_a_a_2a"){
        vparameters.push_back("3.5320797,4.1213427,0.489,0,0.4176,0.5824");  // 001, ternary metal-carbo-nitride prototype (ICSD #95265)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B3C2_hR11_160_6a_3a_2a"){
        vparameters.push_back("3.2481558,12.319746,0.691,0.085,0.8069,0.1911,0.9137,0.3092,0,0.1351,0.8641,0.7406,0.2597");  // 001, ternary metal-carbo-nitride prototype (ICSD #14399)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B3C4_hR15_160_8a_3a_4a"){
        vparameters.push_back("3.2096311,17.15777,0.6852,0.0633,0.7752,0.1494,0.8522,0.2268,0.9365,0.315,0,0.8112,0.1898,0.7212,0.0993,0.9005,0.2788");  // 001, ternary metal-carbo-nitride prototype (ICSD #14401)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B3C2_hR11_166_3c_ac_c"){
        vparameters.push_back("3.2480036,12.324532,0.3091,0.9157,0.1921,0.8654,0.2598");  // 001, ternary metal-carbo-nitride prototype (ICSD #41260)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B3C4_hR15_166_4c_ac_2c"){
        vparameters.push_back("3.2109784,17.150497,0.3149,0.9366,0.224,0.8498,0.1883,0.2788,0.9008");  // 001, ternary metal-carbo-nitride prototype (ICSD #41261)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_cP10_215_e_ab_e"){
        vparameters.push_back("6.32,0.2,0.3");  // 001, ternary metal-carbo-nitride prototype (ICSD #20748)
      }
    }
    // ---------------------------------------------------------------------------
    // Part 3 //DX20210427
    // ---------------------------------------------------------------------------
    if(library=="all" || library == "" || library=="part3"){
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B11CD8E_aP26_1_5a_11a_a_8a_a"){
        vparameters.push_back("4.8558,1.20466246551,1.81123604761,76.515,81.528,71.392,0.1957,0.1206,0.62783,0.2966,0.2336,0.7408,0.0562,0.3094,0.86872,0.1694,0.4162,-0.0197,-0.0682,0.5331,0.0947,0.158,-0.022,0.679,0.45,0.117,0.794,-0.096,0.433,0.816,0.246,0.548,-0.081,0.764,0.632,0.037,0.867,0.409,0.171,0.893,0.248,0.503,0.552,0.42,0.642,0.781,0.158,-0.018,0.345,0.159,0.144,0.035,0.799,0.113,0.14247,0.62758,0.42059,-0.0431,0.2976,0.55444,0.5597,0.2447,0.39414,0.708,0.8369,0.55456,0.3305,-0.052,0.37321,0.3613,0.4537,0.66271,-0.0376,0.1035,-0.05381,0.4056,0.2355,0.05675,0.0303,0.6739,0.17348,0.47392,0.02789,0.47548");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C5_aP18_2_2i_2i_5i"){
        vparameters.push_back("3.16,1.87974683544,2.83227848101,103.9,91.0,92.0,0.67,0.67,0.34,0.36,0.88,0.17,0.743,0.213,0.36,0.245,0.374,0.1,0.244,0.709,0.054,0.208,0.092,0.19,0.735,0.476,0.249,0.562,0.842,0.288,0.77,0.698,0.489");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC8D3_aP26_2_i_i_8i_3i"){
        vparameters.push_back("7.1576,1.03716189784,1.07544707723,115.11511,107.37724,100.55864,0.20773,0.84122,0.1767,0.14441,0.27518,0.25422,-0.03356,0.87658,0.13472,0.72111,0.40891,0.41381,0.80981,0.2953,0.07942,0.74076,0.0321,0.32928,0.26814,0.70828,0.31212,0.22592,0.32857,0.71731,0.38942,0.10112,0.31914,0.43708,0.31682,0.05324,0.23701,0.18228,0.8249,0.68677,0.41955,0.19963,0.6404,0.203,0.4392");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C3_aP28_2_2i_6i_6i"){
        vparameters.push_back("6.56,1.07317073171,1.07317073171,120.0,92.5,101.16667,0.25,0.653,0.431,0.25,0.31,0.764,0.25,0.319,0.431,0.25,0.319,0.097,0.25,-0.01,0.431,0.25,0.653,0.097,0.25,-0.01,0.764,0.25,0.65,0.764,0.25,0.431,0.319,0.25,0.764,0.319,0.25,0.764,0.653,0.25,0.21,0.542,0.25,0.21,0.875,0.25,0.542,0.875");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_aP32_2_12i_4i"){
        vparameters.push_back("7.309,1.02914215351,1.05048570256,88.81,89.08,89.07,0.0007,0.0386,0.21,0.5038,0.5361,0.2181,0.0076,0.466,0.2884,0.4972,0.9638,0.2878,0.2851,0.2574,0.287,0.2204,0.763,0.2232,0.2186,0.2627,0.7258,0.284,0.7583,0.7679,0.2943,0.0422,-0.0002,0.2971,0.5446,0.4982,0.2096,0.482,-0.0072,0.2088,-0.017,0.5051,0.2566,0.0259,0.285,0.2502,0.528,0.2158,0.2438,0.0313,0.7817,0.2499,0.5338,0.719");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C10D3_aP36_2_ah2i_2i_10i_3i"){
        vparameters.push_back("8.102,1.0144408788,1.05801036781,69.15,62.88,67.23,0.65388,0.03726,0.79895,0.71664,0.85636,0.39005,0.004,0.867,0.704,-0.065,0.823,0.624,0.3023,0.7124,0.3858,0.5468,0.6953,0.5643,0.2128,0.8964,0.7589,0.5375,-0.0553,0.6739,-0.0504,0.7241,0.4515,0.7083,0.518,0.2467,0.5988,0.8268,0.0273,0.8667,0.7852,0.1424,0.1894,0.8592,0.1263,-0.0821,0.8725,0.6893,0.56457,0.23365,0.23877,0.79662,0.66737,0.06492,0.17624,0.66966,0.29811");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB10C9D_aP42_2_ae_10i_9i_i"){
        vparameters.push_back("6.141,1.74824947077,0.974759811106,82.26667,107.43333,102.66667,0.898,0.1412,0.2547,0.7185,0.0126,0.2283,0.301,0.2016,0.0667,0.3341,0.127,0.3188,0.3231,0.3785,0.3406,0.6016,0.3937,0.4256,0.8012,0.4011,0.8847,0.857,0.3845,0.162,0.6033,0.1321,0.6671,0.4108,0.1932,0.6922,0.9072,0.152,0.6734,0.2442,0.3172,0.796,0.8061,0.3724,0.6363,0.0444,0.3022,0.3849,0.8176,0.0737,0.1519,0.2887,0.1177,0.149,0.4654,0.4063,0.2975,0.756,0.4161,0.0191,0.435,0.1263,0.6289,0.0133,0.2871,0.6253");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_aP44_2_4i_14i_4i"){
        vparameters.push_back("6.612,1.00862068966,1.82773744707,85.81,89.38,88.57,-0.0521,0.331,0.11666,0.8845,0.0908,0.35915,0.3705,0.7756,0.36947,0.6657,0.828,0.1059,0.6408,0.4893,0.1258,0.6226,0.1401,0.2096,0.2968,0.2956,0.0948,0.4002,0.4255,0.304,0.5879,0.1717,0.4454,0.224,0.095,0.3785,0.2937,0.422,0.5099,0.2862,0.2179,0.687,-0.0374,0.2281,0.5714,0.0752,0.5789,0.6864,0.2414,-0.0859,-0.001,0.3397,0.7841,0.1903,-0.0041,0.6708,0.0797,0.0023,0.0137,0.1857,0.1539,0.8505,0.1168,0.4862,0.3353,0.1761,0.3781,0.2726,0.4051,0.1457,0.3719,0.6179");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B2CD12_aP54_2_12i_2i_i_12i"){
        vparameters.push_back("7.694,1.54873927736,0.756043670393,102.3,102.4,105.9,-0.098,0.375,0.1583,0.1088,0.418,0.275,-0.0233,0.1283,0.355,0.202,0.135,0.442,0.3033,0.1153,-0.095,0.07,0.04,0.7433,-0.0353,0.2655,0.6483,0.1356,0.3861,0.7517,0.4,0.355,0.3917,0.4704,0.3451,0.1817,0.7167,0.1417,-0.0167,0.7667,0.1167,0.7,0.4207,0.873,0.6648,0.3707,0.6081,0.19,0.09,0.2325,0.0606,0.4696,0.8648,0.8845,0.5323,0.8766,0.5422,0.2611,0.8757,0.5844,0.5386,0.6283,0.2804,0.2587,0.6035,0.3158,0.3106,0.5904,-0.0424,0.0135,0.3655,0.2694,0.1112,0.1448,0.333,0.17,0.1038,0.8533,0.0773,0.3149,0.776,0.3614,0.3496,0.2179,0.8166,0.1262,0.9031");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C5D_mP20_4_2a_2a_5a_a"){
        vparameters.push_back("5.449,0.886768214351,1.49330152322,107.19,-0.02933,0.38995,0.31142,0.00673,0.63149,0.44573,0.3031,0.49638,-0.00741,0.56046,0.48968,0.39524,0.02087,0.07409,0.17122,0.43724,0.11111,0.37857,0.39807,0.12076,0.07634,0.32462,0.69627,0.20975,-0.08635,0.46436,0.40379,0.29229,0.0,0.20805");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B6C_mP20_4_3a_6a_a"){
        vparameters.push_back("5.7275,1.03996508075,1.44883457006,90.568,0.2621,-0.0025,0.7488,0.2723,0.4223,0.0123,0.7532,0.0288,0.4761,0.6319,0.41,0.5173,0.5605,0.1618,0.1938,0.0747,0.2744,0.3323,0.1165,0.0484,0.0152,0.0327,0.3053,0.6942,0.574,0.1963,0.8252,0.2459,0.0,0.2519");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B2C2_mP60_4_22a_4a_4a"){
        vparameters.push_back("7.822,1.59805676298,0.991306571209,91.05,0.224,0.099,0.084,0.796,0.889,0.811,0.258,-0.096,0.129,0.839,0.079,0.842,0.193,0.009,0.442,0.852,-0.005,0.561,0.834,0.017,0.197,0.086,-0.024,0.788,0.517,-0.012,0.57,0.991,0.229,0.991,0.088,0.243,0.331,-0.089,0.3,0.665,0.421,0.263,0.257,0.584,0.267,0.722,0.652,0.239,0.082,0.326,0.287,-0.078,0.41,0.01,0.888,0.497,0.007,0.208,0.279,0.162,0.664,0.687,0.842,0.389,0.301,0.841,0.669,0.64,0.155,0.343,0.3226,0.0046,0.0737,0.8943,-0.0084,0.7508,0.6178,0.2582,0.2607,0.382,0.2606,0.7514,0.3082,-0.0024,0.6351,0.7024,-0.0071,0.3585,0.1913,0.2501,0.1533,0.8074,0.2557,0.8462");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mC10_5_b_a_ac"){
        vparameters.push_back("6.7064,0.763241083144,0.627266491709,109.284,0.6259,0.8722,0.103,0.8957,0.5037,0.7367");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC12_5_2c_c"){
        vparameters.push_back("9.357,0.36147269424,0.8327455381,119.76667,0.0948,0.488,0.3928,0.1399,0.067,0.0257,0.3444,0.5,0.3044");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_mC16_5_c_3c"){
        vparameters.push_back("5.91,1.73205076142,1.04286294416,108.64073,0.54333,0.27833,0.045,0.24,0.27778,0.22,0.74,0.44444,0.22,0.74,0.11111,0.22");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C18D4E2_mC62_5_a_2b2c_9c_2c_c"){
        vparameters.push_back("16.8913,0.33393522109,0.494876060457,93.919,0.67531,-0.07889,0.413,0.32143,-0.25821,0.47267,0.08396,0.18216,0.22561,0.33645,-0.2696,0.2451,0.5762,0.1756,0.4677,0.07301,0.1655,0.4508,0.47637,-0.1144,0.2054,0.4428,-0.5435,0.1571,0.3822,-0.2162,-0.0248,0.71142,0.3504,0.4381,0.67462,-0.0515,0.3216,0.61381,0.3236,0.1953,0.64484,0.2009,0.35178,0.4101,-0.29307,0.14827,0.29009,0.20623,0.15752");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C9D2_mC90_5_ab2c_3c_b13c_3c"){
        vparameters.push_back("17.45812,0.39701869388,0.688934432803,133.36555,0.4553,-0.0447,0.5985,0.1667,0.2724,0.8943,0.3333,0.2724,0.6057,0.5,0.171,0.572,0.333,0.307,0.783,0.833,0.023,0.712,0.5173,0.226,0.8882,0.0173,0.226,0.3882,0.8507,0.4435,0.5522,0.3507,0.4435,0.0522,0.184,0.3306,0.1115,0.684,0.3306,0.6115,0.6556,0.1185,-0.0982,0.1556,0.1185,0.4018,-0.0111,0.3102,0.8065,0.4889,0.3102,0.3065,0.3222,0.0714,0.7583,0.8222,0.0714,0.2583,0.3333,0.4508,0.784,0.0833,0.2752,0.3585,0.25,0.4496,0.25,0.4167,0.2752,0.1415");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_mP12_7_2a_4a"){
        vparameters.push_back("8.76,0.503424657534,1.15867579909,125.2,0.7483,0.4938,0.1243,0.246,-0.0059,0.8731,-0.0649,-0.0138,0.5728,0.5585,0.052,0.6747,0.0575,0.4464,0.4209,0.4363,0.5585,0.3235");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C15D4_mP54_7_6a_2a_15a_4a"){
        vparameters.push_back("9.0112,0.812000621449,1.21762917259,107.72,0.3409,0.2658,0.07662,0.6782,0.50471,0.1731,0.3148,0.0044,0.3329,-0.0022,0.2541,0.27358,-0.0033,0.7531,0.31251,0.6554,0.2291,-0.09514,0.3179,0.718,0.0741,0.6731,0.2144,0.4123,0.8207,0.4832,0.8963,0.2666,0.3166,0.3293,0.2631,0.6918,0.3234,0.0643,-0.0059,0.1862,0.7163,0.8148,0.1874,-0.0723,0.5119,0.1535,0.497,0.4776,-0.022,0.5282,-0.0241,0.5029,0.1732,0.0159,-0.0177,0.0011,0.7477,-0.0003,0.4988,0.0147,0.234,0.456,0.5016,0.2321,0.7307,0.1924,0.191,0.8718,0.0697,-0.0515,0.1138,0.4252,0.0116,0.3664,0.5049,0.3381,0.6271,-0.0087,0.6556,0.0263,-0.0323,0.0347,-0.0333,0.5304,0.0206");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B23_mP124_7_16a_46a"){
        vparameters.push_back("13.39,0.602091112771,1.25616131441,106.02,0.0612,0.2916,0.4188,-0.0652,0.2923,0.0817,0.0649,0.7928,0.417,-0.0675,0.7935,0.0847,0.1842,0.2087,0.2479,0.8176,0.206,0.2543,0.1862,0.707,0.2448,0.8156,0.7049,0.2563,0.3154,0.2899,0.0796,0.6851,0.2954,0.4187,0.3156,0.7936,0.0774,0.6845,0.7997,0.421,0.4472,0.2043,0.4015,0.5537,0.206,0.097,0.447,0.7036,0.4012,0.5546,0.7049,0.0981,0.0779,-0.0011,0.405,-0.0728,0.0012,0.0891,0.0607,0.499,0.4296,-0.0486,0.5048,0.0783,0.1929,0.4912,0.2339,0.7976,0.504,0.2613,0.19,-0.0029,0.2574,0.8216,-0.004,0.2478,0.3214,0.0082,0.067,0.6813,0.006,0.4303,0.3152,0.5,0.0882,0.6775,0.5067,0.4152,0.4436,0.4978,0.4074,0.546,0.4885,0.0995,0.4521,-0.0109,0.4152,0.558,-0.0026,0.0935,-0.0002,0.2505,-0.0015,0.0024,0.79,-0.0015,0.0606,0.2467,0.1601,-0.0672,0.2429,0.3282,0.0601,0.7542,0.1589,-0.0673,0.7428,0.3317,0.1265,0.2707,0.3274,0.8668,0.2776,0.1659,0.1341,0.7186,0.3288,0.8756,0.713,0.1692,0.2011,0.2464,0.4929,0.80065,0.2582,-0.001,0.1943,0.7351,0.4921,0.7994,0.7476,0.0018,0.2593,0.2437,0.1578,0.7417,0.229,0.3336,0.2579,0.7718,0.1605,0.7351,0.7823,0.33,0.3226,0.2545,0.3243,0.6708,0.2316,0.1706,0.3225,0.7445,0.3215,0.6702,0.7512,0.1657,0.4085,0.258,0.0033,0.5842,0.262,0.493,0.4024,0.7451,-0.0013,0.5864,0.7398,0.4904,0.4542,0.2477,0.1492,0.5416,0.2337,0.3391,0.4514,0.7444,0.1512,0.5319,0.7539,0.3435");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_mC8_8_a_a_b"){
        vparameters.push_back("7.31,0.682626538988,0.608755129959,114.8333,0.0,0.0,0.486,0.5,-0.083,0.306,0.444");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B4_mC102_8_17a11b_8a2b"){
        vparameters.push_back("15.183,0.534940393862,0.812751103208,107.9,0.055,0.174,0.329,0.28,0.767,0.469,-0.077,0.424,0.754,0.031,0.522,0.169,0.502,0.499,0.702,0.23,-0.095,0.212,-0.055,0.825,0.67,0.72,0.232,0.532,0.077,0.576,0.245,-0.034,0.477,0.831,0.297,0.769,0.095,0.788,0.086,0.382,0.604,0.373,-0.091,0.012,0.597,0.018,-0.087,0.617,0.395,0.626,0.09,-0.013,0.402,-0.018,0.185,0.217,0.111,0.364,0.21,0.111,0.176,0.219,0.334,0.491,0.224,0.331,0.367,0.21,0.477,0.498,0.248,0.0,0.314,0.282,0.888,0.124,0.271,0.888,0.322,0.28,0.666,0.008,0.279,0.669,0.131,0.289,0.532,0.318,0.292,0.277,0.182,0.229,0.722");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5C4D2_mC56_9_3a_5a_4a_2a"){
        vparameters.push_back("5.34,1.73052434457,2.75074906367,93.66,0.8852,0.198,0.2303,0.3747,0.3771,0.2327,0.8818,0.5235,0.2204,0.0726,0.1978,0.0,0.2387,0.475,0.0023,0.7395,0.4214,0.0063,0.0438,0.3596,0.156,0.532,0.537,0.155,0.7299,0.363,0.2878,0.2395,0.197,0.2791,0.5334,0.2053,0.1807,0.2025,0.5416,0.296,0.0257,0.3641,0.041,0.5208,0.5335,0.04");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C9D2_mC68_9_2a_4a_9a_2a"){
        vparameters.push_back("8.91,0.577328843996,1.63782267116,100.5,0.8275,0.4297,-0.0002,0.1628,0.4245,0.0004,0.426,0.318,0.096,0.599,0.41,0.378,0.3,0.346,0.38,-0.001,0.386,0.38,0.1349,0.0011,0.2215,0.6809,0.0027,0.2362,-0.0881,0.3075,0.2371,0.0214,0.3111,0.0788,0.7143,0.238,0.0794,0.3271,0.2379,0.0748,0.164,0.1153,-0.066,0.7853,0.1207,-0.065,0.4772,0.0432,-0.0651,0.0632,0.2796,0.1897,0.7342,0.2656,0.1907");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B36C11_mC212_9_6a_36a_11a"){
        vparameters.push_back("12.577,0.577323686094,3.00548620498,102.81,0.831,0.0001,0.7206,0.2029,0.0017,0.7786,0.2482,0.4974,0.0987,0.2817,0.0025,0.8992,0.4748,0.4942,-0.0638,0.058,0.0037,0.0625,0.528,0.01,-0.044,0.082,0.48,0.667,0.731,0.007,0.892,0.749,0.02,0.791,0.814,0.009,0.007,0.207,0.324,-0.048,0.701,0.203,-0.048,0.053,0.311,0.894,0.561,0.202,0.888,0.152,0.318,0.832,0.618,0.168,0.836,0.372,0.326,0.852,0.859,0.191,0.847,0.432,0.303,0.791,-0.047,0.202,0.784,0.64,0.2,0.012,0.127,0.314,0.006,0.414,0.164,-0.014,-0.093,0.313,-0.013,0.411,0.187,0.667,0.896,0.315,0.668,0.185,0.174,0.652,0.689,0.312,0.646,0.091,0.215,0.713,0.579,0.307,0.707,0.804,0.323,0.044,0.317,0.172,0.048,0.487,0.188,0.611,-0.013,0.33,0.607,0.219,0.004,-0.014,-0.057,0.478,0.832,0.537,0.467,0.846,0.797,0.497,0.712,0.292,0.016,0.609,0.495,0.498,0.649,0.502,0.013,0.05,0.25,0.2342,0.0,0.7808,0.2354,0.0,0.483,0.2476,0.8347,0.5381,0.2682,0.6654,0.0553,0.2539,0.6648,-0.0022,0.2447,0.8349,0.6325,0.0051,-0.0841,-0.0894,0.495,0.0835,0.7428,0.0229,0.8341,0.7913,0.4872,0.6644,0.5166,0.0298,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_mP6_11_e_e_e"){
        vparameters.push_back("4.28,0.848130841121,1.41355140187,112.5,0.219,0.448,0.69,-0.058,0.643,0.312");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mP8_11_3e_e"){
        vparameters.push_back("5.4109,0.692823744664,1.74536583563,97.48,0.762,0.554,0.456,0.174,0.888,0.169,0.285,0.656");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C2_mP18_11_2e_e2f_2e"){
        vparameters.push_back("6.921,0.890044791215,1.08900447912,102.79,0.231,-0.0647,0.6398,0.67278,0.6555,0.0407,0.0291,0.33,0.7014,0.2384,0.0741,0.0515,0.2343,0.635,0.053,0.3147");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2CD6_mP20_11_e_2e_e_2e2f"){
        vparameters.push_back("5.2344,1.54592694483,1.25019104386,106.05,0.1474,0.28824,0.1028,0.7517,0.6149,0.7468,0.6232,0.19855,-0.0057,0.8607,0.6383,0.5644,0.8457,-0.038,0.3089,0.6066,0.4604,0.8474");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C2_mP22_11_3e2f_2e_ab"){
        vparameters.push_back("7.5,1.07466666667,0.669333333333,112.0,0.19,0.31,0.508,0.264,0.88,0.43,0.12,0.588,0.709,0.548,0.19,0.09,0.8,0.688,0.09,0.736");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC6_12_i_a"){
        vparameters.push_back("6.9038,0.477925200614,0.988441148353,122.197,0.5048,0.2294");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4_mC14_12_ai_2i"){
        vparameters.push_back("12.214,0.332405436384,0.427296544948,105.0,0.2147,0.3369,0.4286,0.6864,0.1718,0.8123");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_mC16_12_g_2i_i"){
        vparameters.push_back("11.973,0.458531696317,0.593836131295,118.2,0.2493,0.5357,0.2899,0.1783,0.0907,0.1752,0.6314");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_mC20_12_g_i_ij"){
        vparameters.push_back("6.077,1.73177554714,1.11831495804,107.35,0.33258,0.0556,0.1686,0.7593,0.2497,0.2438,0.1612,0.2516");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC2_mC24_12_ij_h_gi"){
        vparameters.push_back("7.37,1.12075983718,0.432835820896,110.3,0.2066,0.3474,0.237,0.631,0.2857,0.3952,0.107,0.118,0.004");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mC24_12_3i_3i"){
        vparameters.push_back("15.98,0.229536921151,0.596307884856,106.0,0.6632,0.1738,-0.0369,0.1761,0.3479,0.4543,0.2613,0.2076,0.3697,0.0838,-0.0661,0.4116");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B4_mC34_12_b6i_2i"){
        vparameters.push_back("17.64,0.239682539683,0.440646258503,115.15,0.587,0.368,0.257,0.613,0.132,0.162,0.79,0.087,-0.086,0.432,0.409,0.194,0.294,0.2915,-0.0081,0.1947");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2D_mC40_12_ad_gh4i_j_bc"){
        vparameters.push_back("7.9174,0.99483416273,0.999823174274,90.378,0.26,0.759,0.266,0.028,-0.031,0.255,0.24,0.454,0.535,0.247,0.248,0.251,0.249");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B12CD6_mC42_12_i_2i2j_a_ij"){
        vparameters.push_back("9.8607,0.720750048171,0.615950186092,93.758,0.3176,0.6122,0.2372,0.2583,0.2693,0.0083,0.2019,0.1095,0.0209,0.2997,0.2784,0.8839,0.1984,0.3151,-0.0429,0.2067,0.2233");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5_mC48_12_2i_ac5i2j"){
        vparameters.push_back("9.92,1.66129032258,1.0,120.0,0.66667,0.19792,0.66667,0.38542,0.66667,0.66667,0.0,0.25,0.33333,0.08333,0.16667,0.79167,0.16667,0.29167,0.41667,0.25,0.04167,-0.08333,0.25,0.54167");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C4_mC52_12_i_gi3j_2j"){
        vparameters.push_back("8.595,1.51576497964,0.834787667248,115.94,0.147,0.2866,0.138,0.64,0.2849,0.8302,0.1476,0.2269,0.0351,0.3105,0.2569,0.1789,0.1265,0.4083,0.00991,0.1856,0.22381,0.71075,0.11813,0.34438");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C5D24E8_mC82_12_h_i_agh_2i5j_2j"){
        vparameters.push_back("9.8359,1.83460588253,0.536321028071,104.75,0.1765,0.278,0.0878,0.196,0.764,0.1085,0.7155,0.3377,0.2928,0.1119,0.0857,0.218,0.1187,0.1709,0.7244,0.1351,0.2519,0.2069,0.3466,0.1344,0.1005,0.344,0.1188,0.5891,0.2806,0.0839,0.2972,0.2884,0.1711,0.8047");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B2C10D2E2_mC84_12_acghj_bdi_5j_2i_j"){
        vparameters.push_back("7.8713,2.11151906292,0.718559831286,90.0,0.67511,0.67477,0.39281,0.24815,0.23274,-0.03687,0.23438,0.53428,0.26356,0.41042,0.25122,0.25569,0.16153,0.01527,0.25519,0.16127,0.48391,0.00143,0.08917,0.24702,0.02156,0.24925,0.25002,0.52741,0.10004,0.24944,0.13414,0.16612,0.24902");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B4_mC102_12_dg8i5j_4ij"){
        vparameters.push_back("15.489,0.521860675318,0.805474853122,107.71667,0.244,0.064,0.173,0.322,0.277,0.235,0.539,0.081,0.582,0.231,-0.028,0.48,0.827,0.31,0.769,0.086,0.781,0.086,0.383,0.401,0.624,0.09,-0.011,0.4,-0.015,0.188,0.216,0.111,0.373,0.211,0.107,0.176,0.216,0.334,0.495,0.283,0.329,0.366,0.223,0.479,0.318,0.285,0.277");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2CD15E2_mC112_12_2i3j_j_ad_g4i5j_2i"){
        vparameters.push_back("12.017,0.798368977282,0.827910460181,95.03,0.2267,0.061,0.743,-0.041,0.728,0.1728,0.0594,0.6841,0.0965,0.6686,0.4354,0.0104,0.7806,0.70695,-0.04418,0.21077,0.52021,0.072,0.276,0.621,0.079,0.349,0.726,0.047,0.271,-0.027,0.17022,0.26326,0.24936,0.1588,0.3743,0.8871,0.1519,0.0973,0.5987,0.1723,0.1457,0.548,0.1918,0.0236,0.3754,0.0387,0.3351,0.6481");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C11D6E4_mC112_12_e_gi2j_i5j_2i2j_2j"){
        vparameters.push_back("15.31076,1.20830056771,0.348121190588,107.07121,0.75,0.32,0.72,0.9,0.9,0.75,0.45,0.4,0.4,0.32,0.17,0.22,0.32,0.08,0.72,0.4,0.37,0.25,0.4,0.37,0.75,0.4,0.25,0.4,0.25,0.42,-0.05,0.25,0.33,0.45,0.4,0.17,0.9,0.4,0.08,0.4,0.37,0.42,-0.03,0.37,0.33,0.47");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B23_mP62_13_4g_c11g"){
        vparameters.push_back("13.384,0.303466826061,1.26143156007,106.27,0.06397,0.5858,0.41655,0.18514,0.4131,0.24568,0.31501,0.59,0.07903,0.44659,0.4079,0.40204,0.0658,-0.0012,0.4164,0.1911,-0.0029,0.2458,0.3191,0.0069,0.077,0.4469,-0.0082,0.4074,0.0645,0.4943,0.1643,0.1293,0.493,0.3303,0.1989,0.4945,0.496,0.2612,0.5144,0.1632,0.3261,0.4891,0.3282,0.4113,0.5015,0.0045,0.457,0.4908,0.1532");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C2_mP14_14_e_ae_e"){
        vparameters.push_back("7.16,0.959497206704,0.958100558659,126.16667,0.185,0.62,0.233,0.412,0.365,0.151,0.397,0.115,0.331");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_mP16_14_e_e_e_e"){
        vparameters.push_back("6.2953,1.05871046654,0.882880879386,118.138,0.3115,0.0907,0.1334,0.03201,0.11772,0.28577,0.674,0.128,0.528,0.8807,0.1478,0.5318");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP16_14_2e_2e"){
        vparameters.push_back("12.6,0.347619047619,0.944444444444,119.5,0.084,0.152,0.706,0.462,0.858,0.62,0.179,0.366,0.588,0.297,0.632,0.436");  // 001, (part 3)
        vparameters.push_back("5.79,0.905008635579,1.84801381693,117.4,0.3042,0.9143,0.2992,0.2891,0.1626,0.1011,0.235,0.402,0.329,0.232,0.669,0.045");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_mP20_14_3e_e_e"){
        vparameters.push_back("4.029,3.42144452718,2.16827997022,97.33333,0.2754,0.19875,0.263,0.6782,-0.00745,0.32171,0.8203,0.09875,-0.03369,0.2408,0.04976,0.1575,0.7825,0.17081,0.55692");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2D_mP20_14_a_3e_e_d"){
        vparameters.push_back("5.7009,0.995807679489,1.72670455542,125.30487,0.2103,-0.0065,0.265,0.3193,0.7621,0.0287,0.7308,0.7126,-0.0283,0.2523,-0.01532,0.7508");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C3_mP20_14_a_3e_de"){
        vparameters.push_back("5.4139,1.03459613218,1.74764587451,124.72045,0.8828,-0.0442,0.7812,0.6828,0.6733,0.455,0.2239,0.7313,0.0617,0.7353,0.5508,0.7474");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C4_mP22_14_e_c2e_2e"){
        vparameters.push_back("6.238,0.819365181148,2.170246874,97.217,0.51,0.703,0.115,0.13,0.35,0.185,0.07,0.05,-0.085,0.186,0.225,0.049,0.796,0.113,0.203");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C4D_mP22_14_2e_e_2e_a"){
        vparameters.push_back("4.294,1.78854215184,3.03446669772,87.26667,0.167,0.207,0.055,0.849,0.382,0.39,0.281,0.046,0.341,0.261,0.338,0.082,0.764,0.3,0.327");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2_mP22_14_a4e_e"){
        vparameters.push_back("6.213,1.01239336874,1.37719298246,94.76,0.4044,0.5381,0.2682,0.0889,0.2101,0.293,0.3891,0.3069,-0.0014,0.2159,0.8852,0.0417,0.2646,0.8851,0.3335");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_mP24_14_4e_e_e"){
        vparameters.push_back("4.6916,2.31309148265,1.1908943644,96.04,0.7933,0.0502,0.1358,0.867,0.1822,0.5371,0.534,-0.0106,0.4976,0.3792,0.2017,0.2513,0.6445,0.1136,0.3696,0.8904,0.8495,0.0941");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mP24_14_ab_4e_e"){
        vparameters.push_back("9.0306,0.736208003898,1.41740305185,96.88,0.86659,0.2934,0.04766,0.09492,-0.0538,0.18234,0.5056,0.1517,0.1707,0.6635,0.7336,0.07207,0.2128,0.4492,0.1813");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD3_mP24_14_e_e_e_3e"){
        vparameters.push_back("3.51,2.76638176638,2.29344729345,111.85,0.2098,0.237,-0.0768,0.74205,0.2539,0.1773,0.4274,0.0047,0.7145,0.1896,0.3668,-0.0709,-0.0117,0.1629,0.7946,0.4958,0.1707,0.06");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP24_14_3e_3e"){
        vparameters.push_back("7.02,1.67948717949,0.968660968661,112.0,0.206,0.013,0.005,0.089,0.115,0.095,-0.035,0.041,0.177,0.369,-0.071,0.212,0.27,0.188,0.312,0.81,0.13,0.298");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_mP24_14_e_4e_e"){
        vparameters.push_back("7.83,1.02273307791,1.21392081737,139.39,0.0899,0.8839,0.2331,0.4059,0.5186,0.3818,0.343,0.217,0.3158,0.2414,0.4109,0.0667,0.6985,0.351,0.4181,0.4242,0.3572,0.2869");  // 001, (part 3)
        vparameters.push_back("6.5034,1.08720054126,1.27324630193,126.58575,0.18086,0.16033,0.28154,0.8026,0.0077,0.2503,0.8835,0.3315,0.3799,0.673,0.1071,0.4748,0.4176,0.2168,0.1277,0.6926,0.1639,0.3047");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_mP24_14_e_e_4e"){
        vparameters.push_back("5.64,1.47695035461,1.57941666667,126.99641,0.0792,0.4672,0.8353,0.5965,0.3132,0.3392,0.7908,0.1508,0.4079,0.2708,0.2767,0.2847,0.5599,0.381,0.1555,0.7557,0.441,0.5037");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_mP28_14_abe_4e_e"){
        vparameters.push_back("6.987,0.852010877344,1.53753112924,128.36609,0.5082,0.705,0.7605,0.6549,-0.031,-0.0981,0.408,0.5641,0.0995,-0.051,0.717,0.896,0.261,0.7826,0.8405,0.1649,0.5773,-0.0999");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BCD_mP28_14_4e_e_e_e"){
        vparameters.push_back("4.284,3.34990662932,3.16137721755,102.62529,0.4857,0.7807,0.2292,-0.06,0.0955,0.1866,0.625,-0.0116,0.3824,0.7875,0.8953,0.0398,0.264,0.708,-0.043,0.7052,-0.0592,0.2041,0.1442,0.3409,-0.0791");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C2D8_mP30_14_e_ce_e_4e"){
        vparameters.push_back("4.99995,1.16452364524,2.06746667467,92.2103,0.3308,0.7994,0.3192,0.2508,-0.0033,0.0834,0.182,0.3,0.3709,0.0975,0.8972,0.3318,0.0762,0.3126,0.4451,0.4518,0.7098,0.4183,0.4339,0.7949,0.2065");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_mP32_14_2e_5e_ab"){
        vparameters.push_back("7.9137,0.687529221477,1.44663052681,108.803,0.1592,0.0621,0.32,0.3462,-0.0015,0.6414,0.4447,0.3202,0.3267,0.8802,0.2831,0.3532,0.6197,0.3043,0.072,0.0163,0.261,0.1295,0.2507,0.1873,-0.0047");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_mP32_14_2e_5e_e"){
        vparameters.push_back("9.16,0.770742358079,0.745633187773,107.58,0.02458,0.12451,0.23428,0.38547,0.146,0.08372,0.1163,0.3782,0.4513,0.2968,0.0698,0.3547,0.3683,0.4587,0.248,0.5941,0.2681,0.4507,0.8839,0.3639,0.0059,0.298,0.5876,0.0402");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_mP32_14_6e_2e"){
        vparameters.push_back("7.3271,1.03238661954,1.44716463539,133.21984,0.782,-0.03,0.782,0.781,0.536,0.779,0.005,0.736,0.723,0.472,0.743,0.258,0.2768,-0.028,-0.0002,0.288,0.502,0.0,-0.03,-0.026,0.717,0.4649,-0.033,0.2189");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2D_mP32_14_e_4e_2e_e"){
        vparameters.push_back("8.48,1.41580188679,1.45127830189,129.09091,0.2507,0.0067,0.5026,0.3935,0.8371,0.6395,0.114,0.1778,0.3679,0.2208,-0.0834,0.313,0.2767,0.0954,0.6901,0.5843,0.1353,0.0613,0.0171,0.364,0.5512,0.2671,0.4723,0.5181");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_mP32_14_4e_4e"){
        vparameters.push_back("9.909,0.97436673731,0.858007871632,97.29,0.3187,0.6355,0.0432,0.0819,0.5427,0.3252,0.3698,0.3607,0.3431,0.1455,0.3439,0.1643,0.1645,0.7187,0.1923,0.2537,0.4782,0.5099,0.4703,0.5276,0.2192,0.1964,0.4483,-0.0492");  // 001, (part 3)
        vparameters.push_back("6.56,2.05792682927,1.47941615854,113.75283,0.359,0.524,0.118,0.567,0.36,0.425,0.137,0.373,0.318,0.328,0.339,0.038,0.641,0.508,0.346,0.093,0.524,0.213,0.608,0.275,0.245,0.067,0.285,0.115");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C_mP40_14_7e_2e_e"){
        vparameters.push_back("5.846,2.17122819022,1.4565514882,90.0,0.0463,0.2431,0.1994,0.4573,0.2567,0.2114,0.0448,0.1072,0.388,0.4599,0.1107,0.4101,0.2166,-0.0207,0.2159,0.1875,0.1216,0.0063,0.5572,0.0792,0.1317,0.2397,0.4404,0.1833,0.7626,0.2845,0.4446,0.2718,0.1288,0.2229");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B4C2D6E2FG6_mP54_14_3e_2e_e_3e_e_a_3e"){
        vparameters.push_back("11.33,0.983495145631,0.644571932921,94.76,0.7648,0.8251,-0.0619,0.2402,-0.0953,0.8473,0.2213,-0.0589,0.3316,0.4514,-0.0915,0.1233,0.5826,-0.0741,0.2021,0.4986,0.8473,0.6064,0.6695,0.8204,-0.0303,0.327,0.8561,0.8963,0.2991,-0.0101,0.4229,0.5087,-0.0876,0.2358,-0.0988,0.8216,0.8743,0.1151,-0.0193,0.7552,0.1129,0.8611,0.2058");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B3_mP56_14_11e_3e"){
        vparameters.push_back("17.61,0.523452583759,1.36553094832,100.14,0.3796,0.7206,0.2438,0.092,-0.058,0.2834,0.2204,0.6157,0.3231,0.303,-0.0205,0.3408,0.4356,0.4309,0.3821,0.5254,0.8324,0.4013,0.0663,0.6974,0.4316,0.1617,0.1037,0.4538,0.3158,0.7371,0.4684,0.1517,0.6716,0.1514,0.2388,0.4198,0.6692,0.2234,0.8414,0.2354,0.1795,0.8554,0.3955,0.3835,0.7036,0.3654");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_mP60_14_3e_9e_3e"){
        vparameters.push_back("7.066,1.03594678743,2.18313048401,95.40417,-0.0288,0.6242,0.2482,0.7397,0.3735,0.4011,0.7364,0.8791,0.3987,0.6685,0.6253,0.3,0.3031,0.6241,0.2156,0.0328,0.8603,0.349,0.0348,0.3843,0.3473,0.2388,0.8774,0.5086,0.2347,0.3824,0.5078,0.406,0.8038,0.3642,0.4067,0.4467,0.3633,0.2767,0.1245,0.3906,0.2313,-0.0907,0.4076,0.2313,0.3402,0.4075,0.4432,0.6239,0.3016");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB20C2D14E2_mP78_14_a_10e_e_7e_e"){
        vparameters.push_back("6.3016,1.96450107909,1.4616129237,106.112,0.20966,0.66679,0.06554,0.40162,0.70084,0.22515,0.46879,0.66573,0.07617,0.36184,0.57533,0.16861,0.32899,-0.0962,0.22103,0.1068,0.87132,0.25555,-0.0597,-0.09535,0.73282,0.00217,0.81595,0.85798,0.31563,0.05877,0.89729,0.3095,0.14082,0.02457,0.35953,0.65214,0.13458,0.60057,0.76827,0.41576,0.78118,-0.0758,0.54888,0.63349,-0.07032,0.28057,-0.04314,0.82138,0.39147,0.17726,0.88327,0.17658,0.03043,0.89108,0.83608,0.2821,0.06534,-0.00531,0.74517,0.86102,0.41069");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2CD12E2_mP100_14_8e_2e_ad_12e_2e"){
        vparameters.push_back("9.902,0.962532821652,1.21500706928,95.02,0.387,0.729,-0.07,0.275,0.649,-0.074,0.62,0.735,0.074,0.721,0.648,0.081,0.025,0.72,-0.049,-0.037,0.726,0.043,0.733,-0.008,0.06,0.729,0.014,-0.044,0.29373,0.73232,0.16965,0.25852,0.25782,0.16992,0.1163,0.6276,0.8466,0.8888,0.6261,0.1633,0.5998,-0.0908,0.1493,0.5455,0.1528,0.1791,0.0611,-0.0127,0.1725,0.0973,0.4938,0.1853,0.4342,0.4911,0.1674,0.3741,-0.0227,0.1931,0.3537,0.6684,-0.0342,0.6484,0.6653,0.0431,-0.0131,0.7721,-0.0051,0.7802,0.0103,0.0111,-0.04519,0.49523,0.20732,0.51973,0.00747,0.21252");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_mC4_15_e"){
        vparameters.push_back("2.766,2.91142443962,1.20462762111,92.03333,0.131");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_mC20_15_ad_f_e"){
        vparameters.push_back("8.65,1.00751445087,0.703121387283,130.86424,0.625,0.804,0.445,0.865");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_mC24_15_2f_ce"){
        vparameters.push_back("12.061,0.40096177763,0.446314567615,103.12,0.2851,0.1918,0.0517,0.6746,0.0939,0.4122,-0.0351");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_mC32_15_e_e2f_f"){
        vparameters.push_back("3.496,3.66962242563,4.47397025172,90.53,0.70113,0.32714,0.5053,0.08039,0.13798,-0.00513,0.14565,-0.04913,-0.00793,0.22135,0.11044");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5CD_mC32_15_e_e2f_e_b"){
        vparameters.push_back("6.549,1.32768361582,1.07802717972,113.87,0.8323,0.5714,0.1828,0.1855,0.0663,0.4102,0.1025,0.2893,0.1185");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4D2_mC36_15_f_e_2f_f"){
        vparameters.push_back("9.662,0.657213827365,1.147588491,109.4,0.4733,0.2396,0.3216,0.0546,0.1355,0.2748,0.0933,0.3262,0.4768,0.094,0.129,0.8156,0.1297");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C2_mC44_15_e3f_f_f"){
        vparameters.push_back("7.429,1.12262754072,1.35926773455,111.37,0.0612,0.3984,-0.0189,0.362,0.244,0.1541,0.1056,0.1531,0.8353,0.1138,0.2016,0.0049,0.2058,-0.04958,0.3239,0.51955");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C6D_mC48_15_e_2f_3f_e"){
        vparameters.push_back("6.277,2.41851202804,0.90361637725,114.11,0.1705,0.67273,0.742,0.087,0.766,0.756,0.02,-0.077,0.08319,0.27218,0.59103,0.19997,0.38195,-0.08702,0.79177,0.06826,-0.07831");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C4D4E_mC56_15_e_2f_2f_2f_a"){
        vparameters.push_back("12.07,1.12758906379,0.556752278376,107.54,0.37049,0.0641,0.125,0.0557,0.6475,0.4437,0.0905,0.3588,0.1896,0.3431,0.4018,0.0551,0.0562,0.0997,0.2032,0.0914,0.7393,0.4076,0.1439");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC2_mC64_15_5f_f_2f"){
        vparameters.push_back("14.43,0.466597366597,0.721413721414,122.13,0.089,0.002,0.143,0.118,0.287,0.318,0.297,0.429,0.06,0.298,0.157,0.33,0.485,0.102,0.103,0.181,0.093,0.308,0.037,0.257,0.466,0.359,0.122,0.165");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5CD2_mC72_15_f_5f_f_2f"){
        vparameters.push_back("5.14,1.73151750973,3.60894941634,99.91667,0.0,0.33333,0.0,0.20278,0.5,0.05833,0.20278,0.16667,0.05833,0.025,0.08333,0.17639,0.525,0.08333,0.17639,0.275,0.33333,0.17639,0.20278,0.83333,0.05833,0.76111,0.0,0.14306,0.26111,0.16667,0.14306");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC10D2E4_mC76_15_f_e_5f_f_2f"){
        vparameters.push_back("5.1988,1.73628529661,3.86739247519,95.782,0.0992,0.2506,0.0838,0.0002,0.3872,0.2525,0.0543,0.0366,0.4431,0.4459,0.4178,0.0931,0.1685,0.2475,0.3712,0.1685,0.2509,0.3132,0.3424,0.0422,0.0622,0.4492,0.451,0.2587,0.1355,0.0354,0.4298,0.3646");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD12E3_mC76_15_f_e_b_6f_ef"){
        vparameters.push_back("12.004,1.04406864379,0.533488837054,114.4,0.2599,0.7145,0.2812,0.6525,0.3713,0.4533,0.7152,0.5342,0.0988,0.6375,0.2401,0.3272,0.6633,0.102,0.1213,0.3974,0.3119,0.2251,0.822,0.3172,0.3102,0.5021,0.3735,0.2424,0.8911,0.1325");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B5C4D2_mC112_15_a3ef_5f_4f_2f"){
        vparameters.push_back("5.305,1.73213949105,5.37229029218,97.147,0.16667,-0.16667,0.5,0.0,0.33333,0.0,-0.308,-0.33333,0.039,-0.308,0.0,0.039,-0.006,0.08333,0.114,-0.006,-0.417,0.114,-0.256,-0.16667,0.114,-0.308,0.33333,0.039,0.142,0.0,0.211,0.142,0.33333,0.211,0.142,-0.33333,0.211,-0.269,0.0,0.094,-0.269,-0.33333,0.094");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C2D17E6_mC124_15_f_2f_f_e8f_3f"){
        vparameters.push_back("12.6188,0.584691095825,1.10898025169,103.762,0.17321,0.16817,0.32317,0.4971,0.084,0.4645,0.7449,0.4773,0.0419,0.7486,0.02899,0.34078,-0.07147,0.0727,0.03247,0.62614,0.24036,0.1803,0.1154,0.24045,0.16471,0.57054,0.23677,0.43859,0.24484,0.1548,0.49864,0.05623,0.12854,0.05189,0.87541,0.05237,0.25694,0.43729,0.01209,0.43098,0.74466,0.20415,0.01965,0.64127,0.25166,0.1025,0.86648,0.02964,0.09397,0.36195");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C9D3E_mC144_15_2f_bcdef_9f_3f_ae"){
        vparameters.push_back("23.8927,0.311057352246,0.843270120162,147.42281,0.4941,-0.0074,0.8315,0.4922,0.0362,0.3371,0.0092,0.8816,0.2531,0.2432,0.255,0.8661,0.0088,-0.0767,0.3666,0.49,0.0586,0.0493,0.4955,0.1747,0.0711,0.1973,0.1285,0.0641,0.1942,0.2558,0.2265,0.3221,0.3544,0.562,0.298,0.1203,0.569,0.2938,0.2625,0.7239,0.1738,0.3455,0.5976,0.2004,0.2221,0.0995,0.2923,0.2252,0.1972,0.0015,0.5723");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B16C20D3_mC168_15_ef_8f_10f_ef"){
        vparameters.push_back("14.818,0.803279794844,0.638952625186,97.39,-0.05363,0.50828,0.15443,0.40511,0.035,0.083,0.252,0.716,0.095,0.209,0.601,0.251,0.079,0.406,0.271,0.075,0.309,0.195,0.326,0.288,0.273,0.347,0.273,0.094,0.203,0.346,0.06,0.3,0.402,0.0258,0.4358,0.1357,0.4224,0.0799,0.1936,0.3962,0.4709,0.0521,0.1008,0.2269,0.0033,0.0848,0.0936,0.1832,0.2296,0.1152,0.0982,0.0936,0.1928,0.6827,0.2796,0.0746,0.4019,0.2214,0.3584,0.2573,0.0955,0.2506,0.4077,0.12978,0.11705,0.05632");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2CD12E2_mC200_15_8f_2f_ce_2e11f_2f"){
        vparameters.push_back("22.291,0.42837916648,0.539903997129,117.548,0.2496,0.478,0.0215,0.1942,0.4916,0.023,0.1939,0.0004,0.5205,0.3664,0.1105,0.0332,0.3651,0.3925,0.0298,0.4808,0.036,0.1753,0.481,0.4628,0.1753,0.1379,0.2306,0.4589,0.1408,0.258,0.334,0.3794,0.017,0.2991,0.372,0.4906,0.292,0.4444,0.1248,0.036,0.4419,0.3781,0.0326,0.2991,0.3418,0.3983,0.2731,0.0985,0.3462,0.0301,0.2498,0.1074,0.5483,0.257,0.1144,0.2827,0.2554,0.2006,0.1872,0.2697,0.2437,0.3248,0.0844,0.0375,0.3231,0.417,0.0324,0.1102,0.2518,0.3703,0.4779,0.2528,0.0207,0.2598,0.2429,0.2977");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP40_17_abcd_2e_abcd4e"){
        vparameters.push_back("5.504,1.01162790698,2.81976744186,0.481,-0.025,0.019,0.525,0.519,0.025,-0.019,0.475,0.014,0.005,0.375,0.514,0.505,0.375,0.25,0.25,0.375,0.75,0.25,0.375,0.25,0.75,0.375,0.75,0.75,0.375");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP12_18_2c_c"){
        vparameters.push_back("4.898,1.75091874234,0.888321763985,0.759,0.281,0.173,0.855,0.036,0.727,-0.0304,0.1016,0.1358");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB12C5D2_oP40_18_a_6c_b2c_c"){
        vparameters.push_back("5.7715,1.77352508014,1.13412457767,-0.0434,0.395,0.35,0.248,0.67,0.5,0.354,0.67,0.0,0.259,-0.09,0.09,0.36,-0.04,0.41,0.118,0.08,0.354,0.01,0.65,0.449,0.3291,0.553,0.007,0.3377,0.853,0.1009,0.3376,0.5374");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP16_19_2a_2a"){
        vparameters.push_back("6.038,0.934580987082,1.67969526333,0.4174,-0.0911,0.0318,0.1338,0.6367,0.3313,0.3086,0.1404,0.2838,0.4287,0.402,0.1341");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C_oP20_19_2a_2a_a"){
        vparameters.push_back("4.905,1.04852191641,1.72742099898,0.501,0.36,0.857,0.28,0.79,0.17,0.1155,0.1256,0.0795,0.167,0.3172,0.7198,0.07096,0.64832,0.6242");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP24_19_a_4a_a"){
        vparameters.push_back("9.886,0.623811450536,1.0441027716,0.51429,0.20707,-0.03774,0.6855,0.10957,0.85152,0.33503,0.07337,0.87726,0.52257,0.07458,0.15367,0.50874,0.55281,-0.03214,0.71343,0.31122,0.37466");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC4D_oP48_19_6a_a_4a_a"){
        vparameters.push_back("7.503,1.00119952019,0.99800079968,0.012,0.574,0.587,0.0,0.388,0.547,-0.068,0.504,0.421,0.113,0.487,0.453,0.789,0.348,0.122,0.154,0.298,0.892,0.011,0.4837,0.5016,0.082,0.6443,0.1192,-0.0853,0.3498,0.1145,0.8628,0.5773,0.872,0.1571,0.41,0.8974,0.0,0.5019,0.0004");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP80_19_5a_10a_5a"){
        vparameters.push_back("6.535,2.97888293803,0.74078041316,0.493,-0.0828,0.688,0.431,0.8469,0.7985,0.251,0.8198,0.631,0.0695,0.8712,0.649,0.15,-0.0586,0.553,0.628,-0.058,0.79,0.376,0.847,0.0,0.303,0.811,0.422,-0.043,0.843,0.507,0.024,-0.028,0.564,0.22,-0.061,0.351,0.6,-0.046,0.38,0.539,0.767,0.892,0.173,0.737,0.567,0.875,0.852,-0.056,0.5625,-0.0901,0.42,0.6075,0.8019,0.771,0.191,0.7545,0.737,-0.004,0.8781,-0.074,0.325,-0.0368,0.711");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A14BC11D_oP108_19_14a_a_11a_a"){
        vparameters.push_back("6.706,1.75902177155,1.78183716075,0.109,0.7344,0.7161,0.875,0.74,0.6932,0.298,-0.0819,0.7675,0.251,0.0439,0.7498,0.689,-0.0769,0.7175,0.74,0.0473,0.6778,0.197,0.1238,0.5537,0.007,0.1119,0.48,-0.023,0.8655,0.3955,0.779,-0.0807,0.4437,0.319,0.7843,0.5251,0.396,-0.0866,0.4995,-0.05,-0.0782,-0.0187,-0.002,-0.0288,0.8664,0.04106,-0.07911,0.60442,0.42704,0.18381,0.57465,0.48169,0.35142,0.68668,0.69542,0.18889,0.70501,0.36068,0.17853,0.77191,0.00529,0.7607,0.66624,0.18986,-0.03077,0.74842,0.7752,-0.03247,0.66931,0.07373,0.08217,0.54556,0.8931,0.87229,0.46088,0.29961,0.8655,0.537,-0.0694,-0.01048,-0.0641,0.49094,0.2263,0.68354");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_20_a_bc"){
        vparameters.push_back("8.46,0.563829787234,0.673758865248,0.33333,0.33333,0.16667,-0.16667,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oC24_20_b_2c_a"){
        vparameters.push_back("7.009,1.0,0.999571978884,0.306,0.198,0.179,0.058,0.172,0.433,0.17,0.941");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C2_oC32_20_b_a2bc_c"){
        vparameters.push_back("10.06,0.819085487078,0.741550695825,0.033,0.0,0.23,0.78,0.19,0.0,0.29,0.29,0.2,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_oF88_22_k_bdefghij_k"){
        vparameters.push_back("10.3832,0.999961476231,1.00259072348,0.2,0.2,0.1689,0.0447,0.4295,0.0706,0.63141,0.63137,0.6114,0.1245,0.124502,0.126");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_oP28_29_a_6a"){
        vparameters.push_back("10.632,0.589164785553,0.59471407073,0.03183,0.23877,0.25,0.082,-0.043,0.106,0.188,0.893,0.121,0.284,0.82,0.112,0.962,0.529,0.39,0.883,0.609,0.28,0.807,0.7,0.172");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC30DE20F2_oP220_29_a_a_30a_a_20a_2a"){
        vparameters.push_back("12.57,0.980906921241,0.984884645982,0.006,0.259,0.256,0.047,0.77,0.28,0.21,0.237,0.305,0.199,0.277,0.175,0.196,0.722,0.693,0.193,0.774,0.825,0.076,0.057,0.235,0.048,0.953,0.781,0.07,0.455,0.22,0.059,0.542,0.768,0.026,0.18,0.45,0.024,0.698,0.96,0.025,0.192,0.048,0.021,0.68,0.547,0.006,0.012,0.502,0.136,0.06,0.507,0.247,0.385,0.557,0.243,0.264,0.012,0.225,0.488,0.752,0.183,0.48,0.885,0.004,0.503,0.502,0.139,0.552,0.49,0.249,0.757,0.49,0.244,0.893,0.445,0.249,0.01,0.743,0.202,0.002,0.108,0.12,0.784,0.34,0.033,0.84,0.245,0.083,0.71,0.225,0.119,0.275,0.665,0.017,0.33,0.755,0.08,0.205,0.775,0.045,0.268,0.711,0.136,0.04,0.811,0.24,0.022,0.978,0.064,0.096,0.976,0.199,0.202,0.904,0.169,0.528,0.19,0.225,0.51,0.011,0.061,0.61,0.046,0.231,0.684,0.097,0.153,0.257,0.236,0.145,0.748,0.759,0.012,0.105,0.258,0.003,0.405,0.247,0.007,0.242,0.397,0.0,0.256,0.098,0.228,0.302,0.591,0.154,0.456,0.801,0.235,0.808,0.425,0.174,0.988,0.182,0.057,0.086,0.521,0.053,0.577,0.472,0.163,0.088,0.915,0.171,0.581,0.085");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C_oP24_31_2b_a3b_a"){
        vparameters.push_back("10.711,0.413313416114,0.395387918962,0.728,0.454,0.289,0.0,0.379,0.174,-0.024,0.246,0.671,-0.037,0.359,0.857,0.064,0.221,0.631,0.335,0.365,0.226,0.335");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6CD8_oP34_31_2a_2a2b_a_4a2b"){
        vparameters.push_back("7.76,1.73453608247,0.677835051546,-0.08333,0.5,0.58333,0.0,0.125,0.0,0.375,0.5,0.25,0.75,-0.08333,0.77778,0.58333,0.27778,0.81111,0.40833,0.68889,-0.09167,0.8125,0.1875,0.5,0.8125,0.3125,0.0,0.84167,-0.03056,0.40833,0.84167,0.53056,-0.09167");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B4_oP102_31_17a11b_8a2b"){
        vparameters.push_back("8.158,1.51287080167,1.77151262564,-0.0055,0.8151,0.0914,-0.0431,-0.0947,0.6226,-0.0915,0.4431,0.1394,0.6728,0.8113,0.1564,0.8092,0.836,0.252,0.414,0.407,0.852,0.401,0.1491,0.6837,0.6749,0.59,0.396,0.5271,0.0054,0.5986,0.2149,0.7111,-0.0332,0.1047,0.238,0.4179,0.5881,0.8978,0.0,0.0901,0.5114,0.1977,0.825,0.772,0.3167,0.5986,0.8248,0.4119,0.3137,0.7309,0.518,0.2877,0.0076,0.2135,0.2138,0.0996,0.2138,-0.098,0.2831,0.2401,0.0829,0.4088,0.2534,0.255,-0.0828,0.2251,-0.019,0.0935,0.2214,0.2229,0.55,0.2259,0.5885,0.5422,0.2244,0.2946,0.7352,0.2291,0.422,0.4252,0.2133,0.5174,0.7338,0.2319,0.2763,0.2849,0.2257,-0.0901,0.7346,0.2181,0.5969,0.0984");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A17B47_oP128_32_a8c_a23c"){
        vparameters.push_back("21.61,0.908375751967,0.182832022212,0.579,0.033,0.02398,0.26126,0.586,0.12899,0.11845,0.4177,0.13648,0.39897,0.4329,0.2426,0.25684,0.579,0.28825,0.06556,0.5848,0.38285,0.19333,0.4334,0.38519,0.36637,0.4205,0.46561,0.05514,0.5551,0.0225,0.263,0.027,0.131,0.1206,-0.052,0.1349,0.4038,-0.011,0.2465,0.2532,0.007,0.2882,0.0716,0.003,0.381,0.1938,-0.008,0.3919,0.3568,-0.047,0.4649,0.0579,-0.001,0.0767,0.0429,0.505,0.0612,0.179,0.514,0.093,0.3224,0.521,0.0585,0.4581,0.483,0.2048,0.0683,0.486,0.1788,0.1999,0.503,0.2064,0.3388,0.503,0.1892,0.4792,0.538,0.2969,0.1672,0.49,0.3297,0.2778,0.481,0.3266,0.4217,0.465,0.3764,0.0925,0.494,0.4697,0.1567,0.477,0.4396,0.2746,0.483,0.4595,0.4132,0.481");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_33_a_a_3a"){
        vparameters.push_back("9.422,0.622054765443,0.562619401401,0.3179,0.0751,0.0,0.0,0.0,0.0706,0.1298,0.0556,0.8087,0.4117,0.8806,0.8382,0.3663,0.3439,0.8147");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_oP44_33_2a_7a_2a"){
        vparameters.push_back("13.87,0.365753424658,0.600576784427,0.12551,0.3373,-0.00169,0.12564,0.33739,0.51409,0.2715,0.4769,0.0876,0.2658,0.4857,0.413,0.3457,0.0706,0.2465,0.4211,0.5557,0.2448,0.5472,0.7858,0.0866,0.5456,0.7882,0.4206,0.5988,0.3526,0.256,0.3205,0.3744,0.2505,0.539,0.6253,0.2498");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BCD6_oP48_33_4a_a_a_6a"){
        vparameters.push_back("7.9241,1.43859870522,0.840953546775,0.2971,0.5446,0.6555,0.2977,0.5006,0.0044,0.2304,0.2003,0.3942,0.2327,0.4189,0.348,0.49584,0.22117,0.852,0.0566,0.4371,0.3269,0.306,0.4658,0.8049,0.2888,0.6655,0.701,0.3,0.5194,0.4588,0.2626,0.3109,0.4526,0.3038,0.4165,0.1445,0.2784,0.6167,0.0497");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oP84_33_6a_3a_12a"){
        vparameters.push_back("11.38,0.554217926186,0.993321616872,0.1401,0.2529,0.539,0.251,0.4193,0.0515,0.2131,0.1191,0.7403,0.3222,0.1006,0.9442,0.4658,0.0183,0.48,0.4748,0.3758,0.2477,0.07351,0.05043,0.0,0.03257,0.65995,0.23386,0.25309,0.11288,0.27955,0.2657,0.1944,0.0459,0.2171,0.4704,0.1734,0.0717,0.0297,0.2009,0.4547,0.1715,0.2567,0.2009,0.375,0.4486,0.1456,0.0322,0.5116,0.0061,0.3009,0.5305,0.1734,0.2946,0.6635,0.1069,0.022,0.7981,0.287,0.2046,0.8336,0.4563,0.137,0.9458,0.1504,0.4768,0.9668");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oP10_34_c_a_c"){
        vparameters.push_back("4.71143,1.11290415012,0.689459463475,0.5,0.425,0.0287,0.543,0.2511,0.1533,0.504");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_oC20_36_b_a_b"){
        vparameters.push_back("8.843,0.618907610539,0.546760149271,0.214,0.23,0.218,0.121,0.642,0.1763,0.1509,0.2898");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC5_oC32_36_b_a_a2b"){
        vparameters.push_back("15.69,0.350031867431,0.343084767368,0.68572,0.20003,0.5803,0.4908,0.16758,0.21738,0.25,0.0949,0.1368,0.6482,0.2541,0.4777,0.4621");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B7C2D2_oC60_36_2b_a3b_2a_b"){
        vparameters.push_back("8.7135,1.75222356114,0.524278418546,0.5852,0.5951,0.7553,0.0886,0.0872,0.0993,0.1735,0.0527,0.1628,0.3264,0.2203,0.1517,0.2897,0.1244,0.0,0.2101,0.0431,0.5074,0.2934,0.2093,0.5024,0.3254,0.1144,0.6523");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oC80_36_4a4b_2a3b"){
        vparameters.push_back("12.22901,0.883554760361,0.566194646991,0.0,0.0,0.2345,0.024,0.233,0.4,0.3814,0.714,0.157,0.712,0.409,0.218,0.1732,0.1177,0.518,0.1723,0.1189,-0.1,0.1972,0.2467,0.217,0.1824,0.4975,0.225,0.12,0.059,0.214,0.152,0.344,0.506,0.151,0.343,-0.08");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oC80_36_2ab_2ab_2a5b"){
        vparameters.push_back("10.825,1.69524249423,0.594457274827,0.41703,0.2604,-0.08321,0.2614,0.2483,0.4215,0.7469,0.4208,0.1807,0.414,0.6788,0.4138,0.2498,0.1667,0.25104,0.2497,0.4979,0.0911,0.0999,0.283,0.4206,0.3987,0.2809,0.4248,0.6506,0.0313,0.0984,0.8479,0.0321,0.0896,0.2478,0.4296,0.0937");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6D15_oC46_38_b_b_2a2d_2ab4d2e"){
        vparameters.push_back("3.949,2.58090655862,3.7277791846,0.0,0.2142,0.856,0.519,0.0,0.404,0.226,0.3186,0.062,0.3153,0.3154,0.119,0.105,0.133,0.311,0.305,0.447,0.351,0.178,0.324,0.05,0.323,0.307");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_oC88_40_abc_2b6c_a3b"){
        vparameters.push_back("11.8887,1.07920125834,0.861860422081,0.56117,0.0,0.47105,0.28712,0.3413,0.244,0.4846,0.4559,0.48386,0.7361,0.21144,0.4933,0.2554,-0.0226,0.01169,0.25276,0.22891,0.1103,0.5171,-0.0426,0.0319,0.3571,0.1036,0.1236,0.5306,0.228,0.8759,0.205,0.2201,0.1013,0.1507,0.1898,0.0325,0.2957,0.3879");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B8CD12_oC104_41_a2b_4b_a_6b"){
        vparameters.push_back("11.062,1.01021605496,0.817302476948,0.4041,0.0,0.189,0.0943,0.3126,0.2042,0.9432,0.4952,0.978,0.247,0.804,0.967,0.111,0.649,0.161,0.159,0.741,0.251,0.316,0.193,0.0672,0.0843,0.3082,0.0832,0.9313,0.4991,0.2591,0.0254,0.4033,0.2448,0.1753,0.2235,0.2816,0.8801,0.5802,0.0151,0.1707,0.7588");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oF24_43_b_a"){
        vparameters.push_back("16.49,0.712553062462,0.410855063675,0.0,0.0749,0.202,0.8118");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_oF40_43_ab_b"){
        vparameters.push_back("9.601,1.44839079263,0.58014790126,0.125,0.065,0.116,0.5,0.068,0.196,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4D_oF64_43_b_a_2b_a"){
        vparameters.push_back("10.53,0.991452991453,0.655270655271,0.512,0.0,0.188,-0.0375,0.1355,0.116,0.0345,0.131,-0.0345,0.116,0.8765");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_oF88_43_b_a3b_b"){
        vparameters.push_back("20.645,0.406054734803,0.312036812788,-0.0908,0.1658,0.3646,0.75,0.2453,0.5622,0.2774,0.1446,0.4368,0.0332,0.1619,0.3477,0.4592,0.19906,0.4054,0.2343");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C2D12E3_oF184_43_b_2b_b_6b_ab"){
        vparameters.push_back("18.296,1.01918452121,0.359914735461,0.0,0.0374,0.0937,0.6152,0.0533,0.1604,0.0616,0.0848,0.1932,0.1731,0.2208,0.0307,0.6179,0.0226,0.0686,0.866,0.0701,0.1818,0.6099,0.0983,0.0351,0.5003,0.2065,0.1529,0.7261,0.1804,0.2273,0.3903,0.0561,0.1896,0.1112,0.1533,0.2113,0.6231");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI8_44_a_a_c"){
        vparameters.push_back("5.384,0.661218424963,1.03324665676,0.88,0.4147,0.8059,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI8_44_a_a_d"){
        vparameters.push_back("3.528,1.74943310658,1.46853741497,0.0,0.4446,0.1701,0.5747");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5CD2_oI40_44_2c_abcde_d_e"){
        vparameters.push_back("8.367,1.28241902713,0.611330225887,0.5912,0.0195,0.374,0.19,0.4256,0.643,0.305,0.041,0.1669,0.1938,0.1465,0.5076,0.1602,0.2055,0.6362,0.2047,0.1613,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oI48_44_6d_ab2cde"){
        vparameters.push_back("7.53776,2.27653838806,1.0,-0.04167,-0.04167,0.25,0.70833,0.75,0.20833,0.09375,0.45833,0.59375,0.45833,0.84375,0.79167,0.34375,0.79167,0.59375,0.125,0.15625,0.125,0.75,0.29167,0.75,0.25,0.04167");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B17C6_oI100_46_ab_b8c_3c"){
        vparameters.push_back("40.92,0.120478983382,0.128787878788,0.0,0.5151,0.5087,0.3564,0.1537,0.1976,0.3485,0.4426,0.2195,0.8082,0.311,0.1589,0.7925,0.3251,0.0956,0.7915,0.2785,0.0362,0.7505,0.3797,0.0827,0.2994,0.4251,0.0249,0.2973,0.1839,0.1419,0.3421,0.1826,0.065,-0.0152,0.5894,0.1278,0.5138,0.5136,0.18965,-0.0151,0.5861");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_oP8_51_e_be_f"){
        vparameters.push_back("5.545,0.730748422002,1.03354373309,0.7423,0.204,0.256");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6D15_oP46_51_f_d_2e2i_aef4i2j"){
        vparameters.push_back("16.635,0.238292756237,0.534295160806,0.0578,0.6874,0.4582,0.6922,0.0524,0.1305,0.3973,0.0593,0.8146,0.1309,0.6149,0.1744,0.8823,0.1636,0.1963,0.0206,0.3386,0.1286,0.3998,0.0583,0.817");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B12CDE6_oP276_52_d4e_18e_ce_de_2d8e"){
        vparameters.push_back("16.119,1.39413114958,0.592530554005,0.2514,0.76061,0.24343,0.3692,0.1162,0.16583,0.57525,-0.01474,0.16928,0.58113,0.48834,0.41462,0.58176,0.25235,0.01987,0.74957,-0.02314,0.604,0.75,0.579,0.412,0.751,0.313,0.737,0.673,0.543,0.25,0.722,0.505,0.734,0.848,0.587,0.303,0.845,0.335,0.514,0.509,0.326,0.417,0.494,0.606,0.768,0.572,0.345,0.264,-0.066,0.682,0.399,-0.01,0.568,0.422,-0.058,0.5,0.359,0.817,0.673,0.429,0.811,0.623,0.565,-0.099,0.679,0.573,0.604,0.695,0.438,0.66,0.561,0.411,0.615,0.475,0.0886,0.65668,0.74995,0.42094,0.58981,0.7477,0.2554,0.70704,0.4384,0.2673,0.66936,0.1506,0.4466,0.51041,0.65,0.2967,0.57117,0.7349,0.4259,0.54662,-0.0646,0.3956,0.66912,0.8444,0.5448,0.60769,0.7631,0.4179,0.63155,0.558");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_oP16_53_eh_ab_g"){
        vparameters.push_back("8.426,0.437930215998,0.970804652267,0.863,0.46,0.869,0.371");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4D2_oP18_53_h_a_i_e"){
        vparameters.push_back("8.0886,0.46309620948,0.916611032812,0.2402,0.3798,0.23998,0.275,0.064,0.099");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C8_oP26_55_h_ag_2g2h"){
        vparameters.push_back("15.074,0.701406395117,0.235305824599,0.21616,0.3812,-0.01109,0.23185,0.15734,0.04391,0.11591,0.21532,0.28335,0.25032,0.11562,0.45742");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BCD2_oP32_55_ghi_f_e_gh"){
        vparameters.push_back("8.258,1.41220634536,1.08077016227,0.22925,0.2312,0.2042,0.0763,0.0819,0.341,0.2489,0.0601,0.1047,0.3047,0.8831,0.186,0.252");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C5_oP32_55_g_fh_eghi"){
        vparameters.push_back("7.36,1.1535326087,0.773097826087,0.25,0.25,0.143,0.172,0.14,0.44,0.09,0.848,0.14,0.44,0.1,0.72,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B11_oP38_55_g3h_a3g2h"){
        vparameters.push_back("11.609,0.97700060298,0.244293220777,0.34,0.21,0.2844,0.3913,0.0429,0.3952,0.1686,0.174,0.13,0.01,0.15,0.32,0.27,0.25,0.4636,0.2962,0.3404,0.0616");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B3C4_oP68_55_2e2fgh2i_adef_2e2f"){
        vparameters.push_back("3.9001,1.0,7.32622240455,0.0695,0.213,0.1402,0.2961,0.4301,0.2871,0.4303,0.3598,0.0699,0.2038,0.2028,0.2971,0.2966,0.2964,0.2721,0.2271,0.1392,0.2266,0.2275,0.3608"); // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oP56_56_2e_e_4e"){
        vparameters.push_back("8.38,1.6491646778,0.597374701671,0.537,0.139,0.624,0.742,0.052,0.365,0.386,0.143,0.123,0.391,0.185,0.633,0.742,-0.009,0.114,0.596,0.112,0.365,0.885,0.112,0.378");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP16_57_a2d_d"){
        vparameters.push_back("3.8,1.94736842105,1.89473684211,0.0,0.78125,0.5625,-0.03125,0.0625,-0.03125");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC_oP24_57_cde_d_a"){
        vparameters.push_back("5.4896,1.4531113378,1.48092757214,0.1726,0.8602,0.03,0.471,0.2013,0.2998,-0.0807,0.081");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP40_57_cd_e_cd2e"){
        vparameters.push_back("5.506,1.01089720305,2.81874318925,0.243,0.696,0.239,0.218,0.191,0.767,0.2566,0.7722,0.6262,0.536,0.532,0.64,-0.034,-0.033,0.61");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_oP8_58_g_g"){
        vparameters.push_back("4.443,2.39522844925,0.886788206167,0.125,0.121,-0.005,0.355");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C6_oP22_58_g_af_gh"){
        vparameters.push_back("4.497,1.20035579275,1.8714698688,0.179,0.56,0.25,0.258,0.316,0.705,0.218,0.139");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oP28_58_4g_3g"){
        vparameters.push_back("15.297,0.804602209584,0.266784336798,0.7111,0.3393,0.8157,0.5236,-0.0325,0.6442,0.4238,0.3974,-0.0967,0.8493,0.7688,0.1386,0.4239,0.156");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC5D2_oP36_58_g_g_3gh_eg"){
        vparameters.push_back("8.306,1.02624608717,0.727546352035,0.24737,0.24952,0.74394,0.2,0.13,0.424,0.6447,0.1079,0.1268,0.104,0.6063,0.13482,0.36423,0.2685,0.3615,0.2778");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C24D8_oP82_58_g_ae2f_2g5h_2h"){
        vparameters.push_back("5.3117,1.76126287253,3.37571775514,0.17748,0.41276,0.24009,0.682,0.225,0.6656,0.1128,0.1502,0.3414,0.1666,0.11279,0.08868,0.6699,0.12137,0.17464,0.1834,0.1211,0.25265,0.433,0.3436,0.11987,-0.0626,0.34888,0.13244,0.17327,0.28464,0.08483,0.67034,0.29399,0.17095");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C3_oP18_59_ef_ab_af"){
        vparameters.push_back("5.745,0.946562228024,0.860226283725,0.5067,0.7629,0.0836,0.8989,-0.0324,0.6045,0.8102,0.43442,0.3832");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_oP28_60_d_c_2d"){
        vparameters.push_back("11.604,0.369269217511,0.535504998276,0.2726,0.1924,0.8296,0.1258,0.0862,0.7268,0.0917,0.2078,0.1485,0.1478");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C6_oP36_60_c_d_3d"){
        vparameters.push_back("14.2661,0.401889794688,0.353950974688,0.3311,0.3389,0.3191,0.2506,0.0963,0.1041,0.0727,0.4189,0.1163,0.099,0.756,0.1236,0.0793");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B12_oP68_60_c2d_6d"){
        vparameters.push_back("12.044,0.681833277981,0.678927266689,0.3939,0.1719,0.8878,0.4768,0.0857,0.7517,0.1119,0.2454,0.7314,0.0906,0.075,0.7736,0.3535,0.1003,-0.0101,0.1032,0.0806,0.508,0.1356,0.0707,0.27,0.3741,0.2538,0.4894,0.3625");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC_oP24_61_c_c_c"){
        vparameters.push_back("6.44,0.944099378882,1.85248447205,-0.006,0.062,0.053,-0.034,-0.083,0.167,0.118,0.23,0.05");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oP28_61_c_2c_a"){
        vparameters.push_back("5.3945,1.03807581796,2.18098062842,0.0042,0.0559,0.3524,0.1961,0.3018,0.0264,0.0673,0.4782,0.3355");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2CD4_oP64_61_c_2c_c_4c"){
        vparameters.push_back("9.762,1.24984634296,0.453800450727,0.10617,0.60704,0.77298,0.00261,0.68871,0.26018,0.23724,0.56757,0.27717,0.3138,0.7228,0.4574,0.0376,0.68766,0.61914,0.1012,0.60302,0.08204,0.18691,0.5345,0.61701,0.33976,0.67302,0.296");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7BCD_oP80_61_7c_c_c_c"){
        vparameters.push_back("12.42,1.02737520129,1.08776167472,0.1062,0.5507,0.5424,0.2364,0.3205,0.5231,0.0098,0.2107,0.4306,0.0092,0.3285,0.6506,0.2852,0.5483,0.3103,0.2769,0.3103,0.2522,0.1207,0.4706,0.1548,0.1075,0.4124,0.3564,0.1875,0.4353,0.2819,0.0629,0.3616,0.4986");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oP80_61_2c_6c_2c"){
        vparameters.push_back("18.2,0.486813186813,0.285714285714,0.13,0.33,0.37,0.13,0.96,0.37,0.06,0.14,0.2,0.06,0.5,0.2,0.05,0.75,0.05,0.19,0.35,0.06,0.19,0.01,0.05,0.2,0.75,0.3,0.03,0.65,0.29,0.22,0.85,0.04");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_oP16_62_c_c_c_c"){
        vparameters.push_back("9.94,0.616700201207,0.855130784708,0.576,0.748,0.192,0.381,0.376,0.561,0.372,0.074");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_oP20_62_c_3c_c"){
        vparameters.push_back("8.96,0.443080357143,1.65959821429,0.665,-0.054,0.784,0.785,0.667,0.504,0.526,0.102,-0.07,0.18");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_oP20_62_c_c_cd"){
        vparameters.push_back("6.4213,0.842804416551,1.425988507,0.7568,0.4166,-0.0848,0.7548,-0.0893,0.8902,-0.0849,0.4492,0.6866");  // 001, (part 3)
        vparameters.push_back("7.65,0.762091503268,0.933333333333,0.145,0.14,0.02,0.68,0.085,0.3,0.175,0.06,0.06");  // 002, (part 3)
        vparameters.push_back("5.7404,0.864295171068,1.38789979792,0.4138,0.2622,0.2597,-0.085,0.4038,0.4225,0.414,0.4736,0.181");  // 003, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_oP24_62_4c_2c"){
        vparameters.push_back("15.985,0.395308101345,0.687832342821,0.00065,0.6743,0.20442,0.40591,0.25609,0.73979,0.43929,0.41364,0.16182,0.06288,0.34361,0.09158");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_oP24_62_a_2cd_c"){
        vparameters.push_back("8.409,0.797835652277,0.574741348555,0.1293,0.7353,0.3646,0.4385,0.18363,0.44979,0.1328,0.0674,0.3083");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oP28_62_2c_2cd_c"){
        vparameters.push_back("7.476,0.770866773676,1.34711075441,0.82623,0.21062,0.51104,0.20406,0.4621,-0.0834,0.2037,0.0582,0.26702,-0.08029,0.1991,0.0412,0.8522");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC_oP28_62_3cd_c_c"){
        vparameters.push_back("7.371,0.850495183829,0.960792294126,0.7187,-0.0104,0.0399,-0.0201,0.3719,0.1647,0.8761,0.8669,0.1658,0.2327,0.1251,0.5733,0.2569");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_oP28_62_c_4c_2c"){
        vparameters.push_back("11.44,0.328671328671,1.23426573427,0.184,0.334,0.695,0.728,0.076,0.184,0.274,0.492,-0.049,0.595,0.355,0.062,0.537,0.614");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_oP28_62_c_6c"){
        vparameters.push_back("13.481,0.228766411987,0.673243824642,0.39692,0.58022,0.81158,0.66343,0.7465,0.59084,0.6767,0.50317,0.52408,0.35102,0.47071,0.24672,0.42026,0.15152");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C_oP32_62_bc_3cd_c"){
        vparameters.push_back("7.675,0.752456026059,0.975335504886,0.1548,0.8583,-0.0911,0.3602,0.0661,0.6434,0.5009,0.4767,0.8404,0.1532,0.7232,0.0145,0.1255");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_oP32_62_ab4c_2c"){
        vparameters.push_back("7.57,0.969749009247,1.02430647292,0.231,0.473,0.222,0.029,0.496,0.238,0.485,0.724,0.471,0.469,0.47,-0.032");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2D_oP32_62_2cd_b_2c_a"){
        vparameters.push_back("12.05,0.755186721992,0.684647302905,-0.05,0.21,0.05,0.79,0.2,0.34,0.825,0.6,0.25,-0.04,0.12");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC2D_oP32_62_2cd_c_d_c"){
        vparameters.push_back("12.05,0.758506224066,0.683817427386,0.308,0.888,0.766,0.099,0.54,0.441,0.512,0.0,0.444,-0.052,0.207,0.685,0.482,0.37");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C4D_oP36_62_d_d_2cd_c"){
        vparameters.push_back("8.7935,0.954762040143,0.528674589185,0.468,0.7957,0.2436,0.4566,0.55945,0.1022,0.63102,0.5824,0.40358,0.75265,0.55741,0.0981,0.4893,0.5923,0.711");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3D3_oP36_62_c_ac_cd_cd"){
        vparameters.push_back("6.03,1.13847429519,1.51243781095,0.8518,0.5556,0.3094,0.2447,0.1951,0.5148,0.3502,0.5018,-0.0669,0.4666,0.7279,-0.0594,0.5651,0.7879");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2_oP36_62_c_2c2d_d"){
        vparameters.push_back("11.068,0.678080954102,0.485905312613,0.042,0.54,0.146,-0.033,0.878,0.838,-0.024,0.035,0.225,0.213,0.049,0.383,0.1412,-0.0056,0.0376");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C4D_oP40_62_d_cd_2cd_c"){
        vparameters.push_back("10.2718,0.851613154462,0.458575906852,-0.0923,-0.0076,0.7204,0.7617,0.574,0.2793,0.7196,0.4195,-0.0318,0.0834,0.7295,0.633,0.4305,-0.011,0.7907,0.1034,0.269");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C3_oP44_62_2c_2c2d_3c"){
        vparameters.push_back("9.794,0.58627731264,1.3977945681,0.1324,0.4084,0.1792,0.736,0.8837,0.1626,0.3229,-0.0759,0.0298,0.1757,0.0907,0.0289,0.3015,0.0293,-0.0807,0.5411,0.779,0.6492,0.541,-0.076");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C2_oP44_62_3c2d_2c_d"){
        vparameters.push_back("13.655,0.596045404614,0.367337971439,0.1524,0.0613,0.0782,0.5588,0.8969,0.3386,0.1821,0.3711,-0.041,0.6214,0.2329,0.0861,0.4875,-0.0477,0.0817,0.7994,0.3745,-0.0096,0.3395");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8BCD_oP44_62_4c2d_c_c_c"){
        vparameters.push_back("16.42,0.490864799026,0.543848964677,0.2585,-0.0805,0.0202,0.2048,0.2249,0.2938,0.4644,0.8699,0.0706,0.8774,0.0742,0.7142,0.145,0.0797,0.1314,0.4442,0.6327,0.1405,0.5383,0.0597");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC8D2_oP52_62_d_c_2c3d_d"){
        vparameters.push_back("8.037,0.960781386089,1.09120318527,0.61445,0.38657,0.4863,0.664,0.8162,0.4282,0.7411,0.4206,0.4192,0.80709,-0.00324,0.06797,0.87368,-0.04233,0.36496,0.60035,0.0782,0.31351,-0.05333,-0.05574,0.1925");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oP56_62_8c_6c"){
        vparameters.push_back("12.428,0.254103636949,1.64467331831,0.1939,0.2262,0.3843,0.8518,0.1218,0.6891,0.298,0.0835,0.4007,0.5269,0.0952,0.4665,0.2977,0.3863,0.0095,0.8272,0.2464,-0.027,0.4824,0.4186,0.3206,0.6371,0.0382,0.5805,0.0155,0.2745,0.1994,0.7941");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B15C5D_oP96_62_cd_3c6d_3cd_c"){
        vparameters.push_back("13.36,0.782934131737,0.504565868263,0.4705,-0.0482,0.771,-0.003,0.527,0.527,0.778,0.337,0.7207,0.3777,0.7041,-0.0519,0.5046,0.4129,0.6038,0.1805,0.1484,-0.0018,0.8394,0.396,0.509,-0.06,0.357,0.52,0.739,0.452,0.528,0.798,0.297,0.696,0.125,0.522,0.807,0.567,0.274,0.809,0.541,0.3975,0.5534,0.8224");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B4C4DE8F2_oP108_62_4c2d_2d_2cd_c_4c2d_d"){
        vparameters.push_back("16.6959,0.695320408004,0.519864158266,0.0619,0.3937,0.1418,0.8446,0.2168,0.3978,0.2602,0.6753,0.36535,0.025,0.46252,0.52249,0.1375,0.59643,0.021,0.2866,0.2614,0.2957,0.326,0.7133,0.6401,0.5221,0.0391,0.1388,0.6648,0.1731,0.0725,0.5689,0.349,0.08,0.367,0.412,0.017,0.347,0.14835,0.04756,0.18788,0.0123,0.5809,0.2893,0.3068,0.0216,0.05,0.3888,0.0593,0.409");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B3_oP112_62_8c4d_4c4d"){
        vparameters.push_back("11.797,0.825548868356,2.22683733152,0.086,0.627,0.041,0.5,-0.034,0.282,0.238,0.331,-0.048,0.825,0.239,0.807,0.192,0.104,-0.032,0.015,-0.057,0.571,0.0148,0.2567,0.107,0.87,0.163,0.018,0.227,0.134,0.5875,0.304,0.634,0.6443,-0.01,0.134,0.7527,-0.039,0.634,0.8663,0.1592,0.072,0.5139,0.8527,0.572,0.6298,0.1723,0.072,0.7612,0.0602,0.572,-0.0645");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD8E3_oP112_62_d_2c_d_4c6d_3d"){
        vparameters.push_back("12.7334,1.07039753719,0.576962947838,0.00647,0.58561,0.0134,0.38486,0.12005,1e-05,0.37761,0.76529,0.37,0.26883,0.05713,0.48974,0.49472,-0.00042,0.34273,0.10285,0.06934,0.50192,0.21704,0.1224,0.21004,0.06068,0.06612,0.00658,0.23777,0.11405,0.84938,0.43489,0.06457,0.78946,0.30497,0.1234,0.52074,0.41552,0.063,0.22624,0.15908,0.1378,0.01543,0.3426,0.13635,0.73115,0.33173,0.13741,0.30525");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5C22D2E8_oP156_62_d_c2d_2c10d_2c_4d"){
        vparameters.push_back("18.544,0.972066436583,0.284836065574,0.1255,0.893,0.7975,0.4595,0.9536,0.778,0.1831,0.5573,0.0697,0.227,0.1239,-0.0094,0.3897,0.1251,0.1634,0.393,0.125,0.0731,0.8196,0.1825,0.1635,0.0599,0.0682,0.1636,0.7275,0.1857,0.0773,0.563,0.0631,0.0772,0.2209,0.1868,-0.0012,0.0743,0.0668,-0.0068,0.7092,0.1973,0.8822,0.3321,0.0504,0.8882,0.0565,0.2003,0.869,0.8288,0.0486,0.8596,0.5504,0.2305,0.8345,0.567,0.0186,0.8336,0.2743,0.2271,-0.0798,0.0637,0.0245,-0.0817,0.7754");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB22C23D2E2_oP200_62_c_11d_3c10d_d_d"){
        vparameters.push_back("14.0135,1.47800335391,0.499225746601,0.6251,-0.0499,0.7409,0.2149,0.4998,0.2167,0.6263,0.5757,0.226,0.783,0.722,0.697,0.133,-0.088,0.721,0.154,0.715,0.533,0.1191,0.835,0.508,0.165,0.701,0.537,0.7851,0.733,0.63,0.2056,0.544,0.772,0.174,0.558,0.8972,0.155,0.473,0.568,0.8819,0.575,0.6483,0.836,0.522,0.6235,-0.0449,0.2482,0.625,0.1278,0.2472,0.6252,0.044,0.5733,0.5373,-0.044,0.7451,0.7118,-0.0443,0.7553,0.6206,0.0452,-0.0774,0.7368,0.1647,0.851,0.5126,0.1645,0.8438,0.8257,0.1614,0.4745,0.5759,0.8379,0.5243,0.6245,0.0003,0.7478,0.625,0.0412,0.2498");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_oC16_63_c_3c"){
        vparameters.push_back("3.78779,4.28999759754,1.0,0.159,-0.0565,0.5565,0.3155");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_oC20_63_c_f_2c"){
        vparameters.push_back("3.07,4.08143322476,1.22475570033,-0.3137,0.2842,0.0724,0.0143,0.3663");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2CD_oC20_63_b_f_c_c"){
        vparameters.push_back("3.863,2.79730779187,1.79911985503,0.7719,0.05275,0.66585,0.5611");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B_oC24_63_c2f_c"){
        vparameters.push_back("3.9875,3.64388714734,3.4407523511,0.316,0.633,-0.067,0.151,0.209,0.434");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A_oC24_63_3f"){
        vparameters.push_back("3.82236,2.79950083195,3.30314256114,0.0284,0.5903,0.2435,0.5551,0.5705,0.3412");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B3_oC32_63_cfg_ce"){
        vparameters.push_back("9.20101,0.781979369656,1.0619486339,0.0254,0.6251,0.2018,0.3147,0.451,0.2219,0.2863");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C6DE2_oC52_63_g_e_fh_c_f"){
        vparameters.push_back("8.275,1.12422960725,1.43697885196,0.33496,0.1862,0.1287,0.5177,0.31221,0.54036,0.22773,0.07717,0.1531,0.3602,0.4531");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC12D2_oC76_63_eg_c_f3gh_g"){
        vparameters.push_back("16.52,0.473365617433,0.340799031477,0.611,0.328,0.778,0.0,-0.092,0.25,0.094,0.0,0.25,0.0,0.406,0.0,0.667,0.375,0.658,0.173,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B3C4_oC68_64_2dfg_ad_2d"){
        vparameters.push_back("28.573,0.13649599272,0.13649599272,-0.0695,0.787,0.8598,0.7039,0.5699,0.7971,0.2972,0.6392,0.2271,0.2279");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C7_oC88_64_ef_df_3f2g"){
        vparameters.push_back("7.17,1.64993026499,2.05020920502,0.25,0.08,0.25,0.08,0.36,0.3,0.08,0.27,0.41,0.13,0.35,0.49,0.21,0.21,0.16,0.21,0.49,0.34");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8_oC22_65_ag_bd2gh"){
        vparameters.push_back("20.66999,0.185437922321,0.189985578126,0.1836,0.094,0.71,0.189");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B8CD2_oC26_65_h_r_a_i"){
        vparameters.push_back("8.18099,1.00314142909,0.458990904524,0.2133,0.2595,0.045,0.312,0.158");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2D_oC24_67_m_a_n_g"){
        vparameters.push_back("11.47,0.659982563208,0.346992153444,0.458,0.608,0.806,0.386,0.348");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_oF56_70_g_h_a"){
        vparameters.push_back("5.8596,2.09987029831,1.67537033245,0.4414,-0.0203,0.0572,0.2137");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_oI16_71_g_i_eh"){
        vparameters.push_back("7.09,1.67136812412,0.764456981664,0.256,0.33,0.148,0.249");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C3_oI20_71_a_in_cj"){
        vparameters.push_back("5.6333,0.99893490494,1.41336339268,0.2192,0.2485,0.2335,0.2166");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B5_oI28_72_j_bfj"){
        vparameters.push_back("13.7,0.511678832117,0.439416058394,0.26,0.378,0.255,0.27,0.345");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4D_oI28_74_a_d_hi_e"){
        vparameters.push_back("5.662,1.0259625574,1.54680324974,0.386,0.0164,0.2748,0.2352,-0.0007");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C2D_oI44_74_h_ij_i_e"){
        vparameters.push_back("7.7077,1.04085524865,1.09664361612,0.88832,0.47954,0.73085,0.206,0.533,0.28237,0.47949,0.274,0.3281,0.416");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C3_tP72_77_8d_ab2c2d_6d"){
        vparameters.push_back("7.98,1.22556390977,0.25,0.25,0.75,0.25,0.076,0.076,0.19,0.076,-0.076,0.31,0.576,0.076,0.19,0.424,0.076,0.31,0.424,0.424,0.31,0.424,0.576,0.19,-0.076,0.424,0.31,-0.076,0.576,0.19,0.25,0.25,0.518,0.25,0.25,0.018,0.12,0.27,0.456,0.33,0.23,0.456,0.25,0.25,0.642,0.12,0.23,-0.044,0.38,0.27,-0.044,0.25,0.25,0.142");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB5C_tP14_85_c_cg_b"){
        vparameters.push_back("6.1768,0.695052454345,0.1975,0.8102,0.3125,-0.0554,0.2994");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B6CD_tP44_85_bcg_3g_ac_e"){
        vparameters.push_back("13.51,0.501110288675,0.58,0.08,0.25,0.0,0.008,0.55,0.05,0.78,0.54,0.85,0.58,0.34,-0.03,0.58");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tP32_86_2g_g_g"){
        vparameters.push_back("10.82,0.549907578558,0.551,0.852,0.86,0.594,0.198,0.69,0.625,0.645,0.14,0.574,0.791,0.17");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C_tP32_86_d_3g_c"){
        vparameters.push_back("8.029,0.983185950928,0.0558,0.2226,-0.0982,0.2792,0.582,-0.0828,0.0895,0.0685,0.2219");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tP40_86_g_g_3g"){
        vparameters.push_back("9.7329,0.63255555898,0.0343,0.7576,0.115,0.426,0.236,0.426,0.836,0.121,-0.05,0.094,0.209,0.094,0.847,0.552,0.172");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B11C2_tP68_86_2g_ab5g_g"){
        vparameters.push_back("12.676,0.441858630483,0.184,0.1182,0.4979,0.1107,0.8064,0.603,0.0301,0.1938,0.7009,0.0013,0.4086,0.6937,0.0415,0.6291,0.7905,0.1484,0.8032,-0.0057,0.4559,0.34,0.1868,0.0247,-0.087,0.0707");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6C2D_tI20_87_a_eh_d_b"){
        vparameters.push_back("5.5571,1.42396213853,0.255,0.289,0.227");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C24D12_tI82_87_a_h_2h2i_hi"){
        vparameters.push_back("12.06,0.627860696517,0.134,0.2113,0.4587,0.3483,0.3066,0.1206,0.3388,0.4104,0.0517,0.35,0.2148,0.2293,0.1289,0.3281,0.3374,0.0851,0.206");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B4C20_tI132_88_a2f_f_5f"){
        vparameters.push_back("14.97999,0.492924227586,0.1366,0.046,0.7009,0.0956,0.218,0.4912,0.0878,0.0544,0.1702,0.0928,0.2127,0.2448,0.0824,0.2227,0.7476,0.1889,0.0739,-0.0926,0.0253,0.0446,0.7555,0.1554,0.1088,0.507");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3D2_tP16_90_c_f_ce_e"){
        vparameters.push_back("8.12,0.541871921182,0.55,0.84,0.11,0.34,0.24");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12BC10D_tP96_92_6b_a_5b_a"){
        vparameters.push_back("6.783,2.69615214507,0.2106,0.70943,0.1097,0.8667,0.0394,0.2215,-0.0746,0.0867,0.5702,0.1559,0.0507,0.5371,0.3532,0.0601,-0.0051,0.4453,0.0744,-0.0101,0.2941,0.1149,0.1727,-0.047,0.0528,0.4705,0.2449,0.0561,0.0658,0.3599,0.085,0.6209,0.6203,0.0658,0.9237,0.6731,0.0003");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_tI44_98_f_bcde_f"){
        vparameters.push_back("7.2313,1.41443447236,0.81048,0.188,0.803,0.5041,-0.0033");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP8_100_b_a_c"){
        vparameters.push_back("6.3,0.587301587302,0.0,0.25,0.65,0.5");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_tP10_100_b_a_bc"){
        vparameters.push_back("5.74,0.862369337979,0.0,0.53,0.75,0.62,0.42");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_tP72_103_abc5d_2d_abc"){
        vparameters.push_back("11.22,0.700534759358,0.377,0.1715,0.223,0.4274,0.391,0.0972,0.1233,0.12,0.13,0.2546,0.255,0.44,0.0919,0.3511,0.146,0.3271,0.4759,0.468,0.4059,0.1432,0.15,0.1033,0.2405,0.0,0.3671,0.2818,0.1031");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_tP7_111_f_a_n"){
        vparameters.push_back("6.34,1.0,0.27,0.225");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4CD2_tP16_113_c_f_a_e"){
        vparameters.push_back("6.3397,0.592457056328,0.37123,0.14562,0.6062,0.111,0.03,0.125");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B12C4D4E_tP46_114_d_3e_e_e_a"){
        vparameters.push_back("8.442,0.757995735608,0.49194,0.082,0.21,0.416,0.183,0.287,0.52,0.063,0.227,0.603,0.1081,0.275,0.5083,0.1301,0.0595,0.1309");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A19B15_tP68_114_bc4e_ac3e"){
        vparameters.push_back("7.5,2.0,0.165,0.3198,0.712,0.118,-0.008,0.384,0.2,0.1333,0.436,0.208,0.212,-0.194,0.375,0.179,0.4021,0.1993,-0.02,0.0982,0.2823,0.1593,0.7157,0.0982,0.1423");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC2_tP4_115_a_c_g"){
        vparameters.push_back("3.72645,1.39381985536,0.19");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C_tI8_119_c_e_a"){
        vparameters.push_back("3.753,3.55502264855,0.35126");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B7C2_tI44_119_i_bdefgh_i"){
        vparameters.push_back("7.2312,1.41410830844,0.1974,0.0683,0.1889,0.1941,0.2471,0.3741,0.2471,0.87294");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C8D_tI72_120_c_2i_2i_b"){
        vparameters.push_back("7.99,1.33767209011,0.22729,0.45696,0.13106,0.1229,0.28418,0.13628,0.12563,0.08408,0.07733,0.13758,0.3956,0.17093");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_tI32_121_g_f2i"){
        vparameters.push_back("9.47,0.484582893347,0.355,0.2851,0.5932,0.25,0.2,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C6D_tI44_121_i_i_ij_c"){
        vparameters.push_back("8.995,0.739188438021,0.2953,0.243,0.11412,0.2783,0.4004,0.212,0.3276,0.1456,0.254");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_tI40_122_e_d_e"){
        vparameters.push_back("9.6922,0.918418934814,0.21203,0.2035,-0.0729,0.1785,0.205,0.0464,0.1639");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4D_tI40_122_e_b_e_a"){
        vparameters.push_back("7.4264,0.933292039211,0.14867,0.22713,0.12266,0.1933,0.08283,0.12675");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_tI48_122_cd_2e"){
        vparameters.push_back("9.5965,1.22841661022,0.16953,0.27647,0.03373,0.24952,0.29644,0.14693,0.116,0.39956");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8BC4D_tI56_122_2e_b_e_a"){
        vparameters.push_back("7.4997,1.00662693174,0.25,0.15,0.125,0.498,0.589,0.063,0.0843,0.1466,0.1151");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC_tP5_123_cg_a_d"){
        vparameters.push_back("4.19,1.89498806683,0.294");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP6_123_d_eh_a"){
        vparameters.push_back("3.61,1.76454293629,0.215");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tP7_123_j_e_a"){
        vparameters.push_back("7.025,0.589893238434,0.2324");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2C_tP11_123_r_f_a"){
        vparameters.push_back("6.41,0.905928237129,0.29,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A8B2C_tP11_123_r_h_a"){
        vparameters.push_back("6.32,0.879746835443,0.1,0.2,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB8C2_tP22_124_a_n_h"){
        vparameters.push_back("6.21,1.77133655395,0.19,0.3,0.11,0.13");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4D2E8_tP32_126_a_b_h_e_k"){
        vparameters.push_back("6.67,1.56371814093,0.53,-0.055,0.06,-0.03,0.34");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B10C2D34E4F9_tP252_126_k_ce2k_f_h8k_k_d2k"){
        vparameters.push_back("15.63,0.756877799104,0.13,0.055,0.11,0.89,0.87,0.81,0.05,0.86,0.17,0.09,0.88,0.83,0.22,0.08,0.87,0.16,0.78,0.22,-0.06,-0.08,-0.07,0.13,-0.02,0.01,0.83,0.82,0.78,0.87,-0.07,-0.05,0.82,0.18,-0.08,-0.1,-0.07,0.06,-0.01,0.83,-0.05,0.19,0.87,-0.09,0.83,0.87");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4D_tP16_127_h_d_i_a"){
        vparameters.push_back("10.302,0.42127742186,0.285,0.194,0.027");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3D2_tP32_127_g_eh_gk_k"){
        vparameters.push_back("8.16,1.08860294118,0.7572,0.3257,0.211,0.3521,0.3726,0.1269,0.1659,0.2594");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B14C5_tP44_128_ac_ehi_bg"){
        vparameters.push_back("7.0138,1.48307622116,0.1711,0.2768,0.0642,0.2477,0.1794,0.5364,0.1198");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC16DE28F8_tP116_128_h_a_2i_b_g3i_i"){
        vparameters.push_back("8.965,1.75872838818,0.1369,0.1094,0.2466,0.4515,0.177,0.0775,0.2362,0.4706,0.1198,0.0846,0.1891,0.2178,0.2636,0.1026,0.0923,0.2131,0.4491,0.0898,0.2256,0.0865,0.19");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_tP8_129_b_c_a_c"){
        vparameters.push_back("4.05,2.23185185185,0.1356,0.6929");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_tP12_129_c_i_a"){
        vparameters.push_back("5.82,0.710652920962,0.48,0.897,0.147");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C6DE_tP26_129_c_j_2ci_a_c"){
        vparameters.push_back("6.98,1.20630372493,0.612,0.343,0.893,0.106,0.584,0.106,0.486,0.392");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A18B10C_tP116_130_2c4g_2c2g_a"){
        vparameters.push_back("8.984,1.28194568121,0.4888,0.2426,0.40482,0.15904,0.8044,0.5759,0.1201,-0.047,0.4874,0.1327,0.8054,0.5766,0.8623,0.8392,0.4596,-0.0398,0.84543,0.47219,0.12549,0.84587,0.47366,0.87655");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C_tP14_136_i_g_b"){
        vparameters.push_back("7.004,0.536122215877,0.2201,0.092,0.3468");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tP20_136_j_dfg"){
        vparameters.push_back("7.63,0.917169069463,0.34,0.2,0.125,0.21");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC4D2E2_tP26_136_fg_a_j_d_e"){
        vparameters.push_back("7.477,1.06125451384,0.2484,0.2161,0.7262,0.0739,0.3178");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB14C2_tP68_136_f_ce2j2k_fg"){
        vparameters.push_back("8.792,1.38648771611,0.116,0.124,0.3572,0.7698,0.0978,0.2942,0.3184,0.255,0.567,0.2245,0.3735,0.1397,0.537,0.1759");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tI10_139_e_ae"){
        vparameters.push_back("3.36301,4.49555309083,0.2,0.4");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C2D_tI18_139_h_d_e_a"){
        vparameters.push_back("6.99,1.25178826896,0.212,0.23");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_tI20_139_ab_eh_d"){
        vparameters.push_back("7.49,1.45126835781,0.288,0.228");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B4C3_tI34_139_c2eg_2e_ae"){
        vparameters.push_back("3.9,7.20512820513,0.068,0.204,0.432,0.296,0.136,0.136");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C3D_tI168_139_egikl2m_ejn_bh2n_acf"){
        vparameters.push_back("15.841,1.1366075374,0.142,0.347,0.142,0.214,0.161,0.673,0.364,0.386,0.114,0.181,0.362,0.16,0.163,0.157,0.276,0.295,0.132,0.293,0.376");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC_tI16_140_h_d_a"){
        vparameters.push_back("5.67,1.20105820106,0.142");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6_tI28_140_a_hk"){
        vparameters.push_back("10.29,0.509232264334,0.4068,0.2141,0.1021");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC2_tI32_140_bl_a_h"){
        vparameters.push_back("8.39,1.70917759237,0.158,0.163,0.363");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC3_tI36_140_cl_b_ah"){
        vparameters.push_back("9.063,1.59439479201,0.66225,0.1421,0.15711");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A31B20_tI204_140_b2gh3m_ac2fh3l"){
        vparameters.push_back("11.076,3.33450794511,0.094,0.1734,0.0756,0.1656,0.1586,0.4035,0.1812,0.0726,0.3417,0.1269,0.1536,0.2109,0.2947,0.4299,0.05,0.2855,0.5774,0.1346,0.2819,0.4125,0.2114");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB11_tI48_141_a_bdi"){
        vparameters.push_back("12.02,0.643926788686,0.123,0.455,0.183");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_tI160_142_deg_3g"){
        vparameters.push_back("12.633,2.01274439959,-0.00073,0.73924,0.24597,0.25789,0.12315,0.13951,0.36959,0.05249,0.11169,0.64224,0.0725,0.11879,0.1061,0.06251");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B2C3D12E4_tI184_142_f_f_be_3g_g"){
        vparameters.push_back("13.723,0.997303796546,0.1242,0.1631,0.3805,0.105,0.3705,0.2186,0.2218,0.1028,0.3633,0.3631,0.2182,0.1051,0.1264,0.1617,0.4118");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C9D_hP28_143_2a_2d_6d_bc"){
        vparameters.push_back("8.84326,0.630483554707,0.36,0.84,0.25,0.77,0.363,0.084,0.229,0.278,0.365,0.737,0.172,0.056,0.839,0.138,0.169,0.418,0.181,0.498,-0.017,0.728,0.198,0.594,0.134,0.502,0.427,0.617,0.143,-0.004");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hP45_144_3a_9a_3a"){
        vparameters.push_back("10.55,0.708056872038,0.4388,0.5668,0.1053,0.0962,0.2042,0.5318,0.742,0.2057,0.1168,0.3371,0.5582,0.0094,0.3964,0.4704,0.228,0.5646,0.6512,0.0664,-0.0011,0.1089,0.6245,0.2309,0.2495,0.5637,0.0581,0.2465,0.3949,0.6922,0.1208,0.2503,0.8721,0.2774,0.0897,0.6541,0.2174,0.014,0.4566,0.5691,0.6236,0.1184,0.2192,0.0,0.7772,0.2214,0.6381");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C_hP12_147_abd_g_d"){
        vparameters.push_back("5.441,1.12718250322,0.67,0.17,0.14,0.4,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD6_hR10_148_c_a_b_f"){
        vparameters.push_back("4.8069,3.32896461337,0.24289,0.49193,0.27933,-0.03947");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B6CD_hR14_148_f_f_b_a"){
        vparameters.push_back("10.59967,1.01326739417,0.94,-0.69,0.14,0.44,-0.19,0.64001");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C6D_hR15_148_f_c_f_a"){
        vparameters.push_back("6.541,1.95887478979,0.28715,0.248,-0.202,0.326,0.296,-0.2014,0.1778");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC8D2_hP12_150_b_a_dg_d"){
        vparameters.push_back("4.7,1.69361702128,0.016,0.222,0.328,0.344,0.317");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B9C3_hP14_150_d_eg_ad"){
        vparameters.push_back("7.37,1.20895522388,0.805,0.32,0.48,0.365,0.2,0.325");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B12C6D_hP21_150_d_2g_ef_a"){
        vparameters.push_back("7.9596,0.518154178602,0.429,0.3114,0.7868,0.4326,0.0988,-0.0926,0.767,0.1113,0.4835");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2CD6_hP30_150_e_c2d_f_3g"){
        vparameters.push_back("8.692,0.707317073171,0.267,0.356,0.817,0.687,0.3586,0.191,0.677,0.355,0.173,0.517,0.822,0.15,0.001,0.255");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hP30_150_ef_3g_c2d"){
        vparameters.push_back("9.82,0.647657841141,0.16,0.59,0.27,0.3,0.62,0.16,0.11,0.23,0.61,0.17,0.34,0.5,0.21,0.8");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2DE3_hR9_155_b_c_c_a_e"){
        vparameters.push_back("4.427,4.2340185227,0.19743,0.27848,0.19063");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hR5_160_a_a_b"){
        vparameters.push_back("6.011,1.35618033605,0.4827,0.0,0.54476,1.11108");  // 001, (part 3)
        vparameters.push_back("5.487,1.66867140514,0.0,0.405,0.565,1.172");  // 002, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2D_hR7_160_a_b_2a_a"){
        vparameters.push_back("5.48,3.87773722628,0.03,0.11,0.54,0.83,0.5,0.0");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B7C_hR11_160_b_a2b_a"){
        vparameters.push_back("8.006,0.857232075943,0.2289,0.0,0.5395,0.14849,0.3895,0.76679,0.833,0.144");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C8_hR13_160_a_ab_2a2b"){
        vparameters.push_back("6.90572,2.42010391386,0.0,0.399,0.636,0.135,0.396,0.814,0.642,0.088,1.138,-0.40801");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B24C_hR28_160_b_2b3c_a"){
        vparameters.push_back("24.817,0.178422049402,0.3146,-0.0551,0.1102,0.5891,-0.00921,0.4515,0.17749,0.8641,-0.2688,0.24779,0.6198,-0.5701,0.08749,0.7609,-0.3733,0.37359");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_hR42_161_2b_b_4b"){
        vparameters.push_back("12.519,1.01629523125,0.5427,-0.4113,0.3732,0.2639,0.3403,0.11859,0.7375,0.01891,-0.62672,0.4582,-0.5851,0.62499,0.6707,0.2584,-0.4056,0.3591,0.2197,0.14411,0.7788,-0.5047,-0.07361");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C_hP9_162_d_k_a"){
        vparameters.push_back("7.906,0.514798886921,0.387,0.0162");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC2_hP9_162_k_a_d"){
        vparameters.push_back("5.295,1.01454202077,0.377,0.2965");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC_hP16_163_i_b_c"){
        vparameters.push_back("5.227,1.90931700784,0.33,0.33,0.15");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6D_hP18_163_d_b_i_c"){
        vparameters.push_back("5.0081,1.92554062419,0.3767,0.031,0.14336");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3_hP4_164_b_ad"){
        vparameters.push_back("3.475,2.45,0.33333");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC2_hP5_164_d_a_d"){
        vparameters.push_back("3.14979,1.51445016969,0.413,0.2203");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_hP9_164_d_a_i"){
        vparameters.push_back("6.73,1.52451708767,0.5,0.135,0.1125");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C_hP9_164_i_d_a"){
        vparameters.push_back("5.784,0.82918395574,0.33,0.139,0.22");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC2_hP9_164_i_a_d"){
        vparameters.push_back("5.62,0.827402135231,0.7,0.148,0.22");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C3_hP12_164_d_ae_i"){
        vparameters.push_back("7.3477,0.720701716183,0.3507,0.8196,0.2492");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_hP12_164_c2d_c2d"){
        vparameters.push_back("4.212,5.44681861349,0.3727,0.1279,0.0419,0.7959,0.2797,0.5596");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A16B2C_hP19_164_2d2i_d_b"){
        vparameters.push_back("2.79596,1.9003812644,0.57451,0.39139,0.83626,0.17157,0.0797,0.1644,0.23653");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C_hR6_166_b_2c_a"){
        vparameters.push_back("3.87818,4.52899555977,0.1118,0.3627");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_hR7_166_a_2c_c"){
        vparameters.push_back("4.036,5.51833498513,0.1453,0.4311,0.2503");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6_hR8_166_a_b_h"){
        vparameters.push_back("5.425,1.8130875576,0.605,1.215");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13B2_hR15_166_b2h_c"){
        vparameters.push_back("6.617,1.82847211727,0.3823,0.8044,1.316,0.9936,1.6714");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C2_hR28_166_2c_2c2h_abh"){
        vparameters.push_back("7.0242,4.95914125452,0.18759,0.35449,0.1232,0.28892,0.19936,0.7283,0.61457,0.13398,1.08544,-0.42793");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B21C24D12_hR62_166_a2c_ehi_fg2h_i"){
        vparameters.push_back("13.80257,1.09220819021,0.2038,0.4065,0.2638,1.1548,0.418,0.75989,0.2515,0.89459,1.0248,-0.67229,0.2024,-0.5087,0.31011,0.1044,-0.1251,0.33381");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hR8_167_e_b"){
        vparameters.push_back("5.198,2.56464024625,0.664");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC6_hR20_167_e_b_f"){
        vparameters.push_back("11.8335,1.00556555541,0.99,0.22,-0.04,0.12");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC3D_hR22_167_f_b_e_a"){
        vparameters.push_back("12.033,1.15208177512,0.62611,0.78748,-0.41951,0.45463");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B3C2_hR28_167_ef_e_c"){
        vparameters.push_back("12.79567,1.4299790476,0.348,0.097,0.58333,0.586,-0.755,0.42601");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_hR42_167_f_ac_2f"){
        vparameters.push_back("7.235,5.41700069109,0.84983,0.35996,-0.25104,0.50497,0.335,-0.429,0.71401,0.4775,0.1905,-0.05849");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A25B21_hR92_167_b2e3f_e3f"){
        vparameters.push_back("25.847,0.339343057221,0.189,0.8054,0.5557,0.3062,0.0702,0.3676,0.4819,0.0723,0.1889,0.235,-0.1163,-0.1187,0.0445,-0.1882,0.2709,0.182,-0.0505,0.0662,0.561,0.1397,0.0436");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_173_b_a_c"){
        vparameters.push_back("5.1815,0.997954260349,0.8907,0.0,0.09358,0.34396,0.1698");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4D_hP14_173_a_b_bc_b"){
        vparameters.push_back("5.146,1.6781966576,0.0,0.3146,0.5358,0.7045,0.3446,0.4031,0.7583");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C7D_hP24_173_a_c_b2c_b"){
        vparameters.push_back("10.31,0.561978661494,0.278,0.024,0.664,0.123,0.357,0.25,0.085,0.25,0.761,0.41,0.526,0.523");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C9D_hP28_173_a_c_3c_b"){
        vparameters.push_back("8.84324,0.630487242232,0.358,0.25,0.3636,0.0854,0.2213,0.17,0.047,0.861,0.192,0.483,0.044,0.142,0.517,0.468");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BCD3E15F3_hP52_173_c_b_b_c_5c_c"){
        vparameters.push_back("12.72,0.407232704403,0.36,0.86,0.26,0.23,0.24,0.5,0.5,0.22,0.05,0.36,0.01,0.36,0.32,-0.03,0.17,0.27,0.26,0.87,0.16,0.24,0.2,0.64,0.36,0.033,0.26,0.26");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_hP20_176_2h_ah"){
        vparameters.push_back("11.064,0.361713665944,0.375,0.861,0.0,0.765,0.717,0.564");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B3C2_hP28_176_hi_af_f"){
        vparameters.push_back("7.17,2.26638772664,0.5718,0.3241,0.4588,0.4472,0.1348,0.3506,0.4074");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B7_hP34_176_c3h_b2h"){
        vparameters.push_back("11.85,0.305738396624,0.8124,0.8157,0.5992,0.8717,0.39,-0.0713,0.7667,0.1058,0.5497,0.1678");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2C9_hP40_176_hi_f_hi"){
        vparameters.push_back("6.436,2.50512740833,0.17175,0.332,0.9098,0.331,0.0895,0.338,0.4243,0.1133,0.3441,0.2817,0.0747");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5BC12D3_hP42_176_fh_a_2hi_h"){
        vparameters.push_back("9.397,0.731935724167,-0.001,-0.00712,0.24227,0.4849,0.3273,0.4667,0.5875,0.36895,0.3985,0.2575,0.3421,0.0705");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD2_hP18_180_f_c_b_i"){
        vparameters.push_back("6.667,1.50802459877,0.3333,0.1521");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC3_hP10_182_c_b_g"){
        vparameters.push_back("5.469,0.942585481807,0.33333");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6_hP18_182_f_b_gh"){
        vparameters.push_back("5.218,1.68282866999,0.054,0.36,0.371");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_hP6_186_b_ab"){
        vparameters.push_back("4.24,3.22405660377,0.375,0.0,0.625");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD_hP8_186_b_b_a_a"){
        vparameters.push_back("3.6648,2.79155752019,0.18,0.0892,0.0079,0.3433");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B8C2_hP26_186_c_ab2c_2b"){
        vparameters.push_back("5.775,1.71688311688,0.3886,0.147,-0.0535,0.5132,0.1461,0.25,0.4861,0.3639,0.1647,0.6354");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC7D_hP26_186_ac_b_a2c_b"){
        vparameters.push_back("5.4317,1.63063129407,0.0629,0.3728,0.6245,0.0,0.1664,0.3126,0.4961,0.3706,0.1616,0.1269");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB6CD7_hP30_186_b_d_a_b2c"){
        vparameters.push_back("7.7192,0.706433309151,0.27671,0.5,0.2381,0.56534,0.08902,0.12232,0.02787,-0.06662,0.26326,0.53382");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B9CD9_hP44_186_c_3c_b_cd"){
        vparameters.push_back("11.73,0.576300085251,0.25,0.13,0.73,0.425,0.49,0.425,0.01,0.215,0.25,0.105,0.53,0.065,0.365,0.75");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B3C2_hP16_187_jk_jk_ck"){
        vparameters.push_back("9.9832,0.423751903197,0.8339,0.0898,0.1676,-0.0873,0.5387");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B_hP24_187_ai2j2kn_j"){
        vparameters.push_back("16.244,0.562977099237,0.20321,0.55085,0.21582,0.74805,0.44928,0.11229,0.81357,0.21992");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C_hP20_190_ac_i_f"){
        vparameters.push_back("6.326,1.82342712615,0.66,0.44,0.33333,0.125");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABCD3_hP36_190_h_g_af_hi"){
        vparameters.push_back("7.1175,1.37153494907,0.449,0.33941,0.032,0.71,0.207,0.891,0.3245,0.3828,0.6354");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B_hP13_191_cdei_a"){
        vparameters.push_back("4.23,1.73286052009,0.29,0.25");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B5_hP18_193_bg_dg"){
        vparameters.push_back("7.861,0.693550438875,0.62,0.29");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B_hP24_193_ack_g"){
        vparameters.push_back("7.12,1.02247191011,0.34,0.33333,0.075");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP10_194_f_bf"){
        vparameters.push_back("4.36079,2.97196838188,0.143,-0.07");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_hP10_194_a_bc_f"){
        vparameters.push_back("2.84,3.80669014085,0.0913");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2_hP10_194_a_f_f"){
        vparameters.push_back("4.0829,4.30943691984,0.82845,0.10706");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B_hP22_194_bhj_c"){
        vparameters.push_back("4.633,0.562702352687,0.115,0.375,0.077");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C2_hP24_194_f_k_bh"){
        vparameters.push_back("8.659,0.787966277861,0.05,0.5231,0.1692,0.0432");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B2C3_hP28_194_hk_f_bf"){
        vparameters.push_back("7.22,2.48337950139,0.837,0.077,0.508,0.824,0.092");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2C9D3E_hP36_194_g_f_hk_h_a"){
        vparameters.push_back("7.39,1.3599458728,0.08,0.47,0.2,0.136,0.125");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_hP60_194_3fk_cdef2k"){
        vparameters.push_back("5.56,4.05575539568,0.14,0.02,0.17,-0.17,-0.05,0.83333,0.1,0.16667,0.05,0.5,0.14");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A12B19C_hP64_194_ab2fk_efh2k_d"){
        vparameters.push_back("5.88,3.91496598639,0.15,0.028,0.19,-0.05,0.182,0.167,0.892,0.167,0.05,0.5,0.15");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4D_cP28_198_a_a_ab_a"){
        vparameters.push_back("7.37,0.258,0.744,0.125,0.0,0.658,0.644,0.0556");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C4D_cP32_198_a_2a_ab_a"){
        vparameters.push_back("7.48,-0.007,0.5,0.75,0.133,0.253,0.556,0.667,0.222");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B4C_cP56_198_ab_2a2b_2a"){
        vparameters.push_back("8.925,0.886,0.0105,0.7621,0.459,0.2113,0.1335,0.1211,0.8719,0.2699,0.4834,0.4706,0.271,0.1892,0.0313");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cI36_199_b_c"){
        vparameters.push_back("7.77,0.4,0.178,0.25,0.403");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B11C3_cP68_201_be_efh_g"){
        vparameters.push_back("9.302,0.38379,0.152,0.59,0.3897,0.599,0.247,0.547");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB3C6D12_cF88_202_a_bc_e_h"){
        vparameters.push_back("10.46,0.195,0.235,0.1");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4BC12_cI34_204_c_a_g"){
        vparameters.push_back("7.832,0.3539,0.1504");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2_cI36_204_d_g"){
        vparameters.push_back("7.6937,0.38587,0.32597,0.1425");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7BC12_cI40_204_bc_a_g"){
        vparameters.push_back("7.312,0.3128,0.1829");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC_cP32_205_d_b_a"){
        vparameters.push_back("8.18,-0.05,0.05,0.225");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C_cP36_205_c_d_a"){
        vparameters.push_back("7.8586,0.34865,0.27772,0.28869,0.4773");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B_cP40_205_cd_c"){
        vparameters.push_back("12.26,0.252,0.125,-0.002,-0.002,0.252");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B2C_cP40_205_bd_c_a"){
        vparameters.push_back("8.2,0.39,0.394,0.218,0.458");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6C6D_cP60_205_c_d_d_a"){
        vparameters.push_back("10.316,0.259,0.195,0.05,0.965,0.19,0.145,0.33");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B6CD6_cP60_205_c_d_a_d"){
        vparameters.push_back("10.96,0.25,0.24,0.0,0.0,0.25,0.25,0.01");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC4_cP84_205_d_ac_2d"){
        vparameters.push_back("9.008,0.37305,0.1189,0.1901,0.3457,0.3336,0.2692,0.1208,0.0906,0.2823,0.0064");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB12CD8E2_cP96_205_a_2d_b_cd_c"){
        vparameters.push_back("12.4,0.239,0.31,0.158,0.014,0.018,0.042,0.136,0.302,0.307,0.224,-0.08");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB24CD20E2_cP192_205_a_4d_b_c3d_c"){
        vparameters.push_back("12.213,0.3343,0.2652,0.585,0.319,0.378,0.486,0.308,0.383,0.555,0.202,0.502,0.588,0.339,0.113,0.2957,0.2783,0.1508,0.0767,0.0403,0.3188,0.1371,-0.0404,0.0573");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB24CD28E2_cP224_205_a_4d_b_2c4d_c"){
        vparameters.push_back("12.135,0.23606,0.3746,0.3061,0.175,0.0365,-0.07275,0.1982,0.0262,0.034,0.104,0.155,0.3,0.0009,0.1821,0.2914,0.26302,0.41979,0.30624,0.2808,0.2,0.3614,0.15173,0.028,-0.01985,0.04692,0.13071,0.30399");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C36D2E20F2_cP252_205_a_c_6d_c_c3d_c"){
        vparameters.push_back("12.314,0.4501,0.5156,0.2341,0.3022,0.3851,0.5072,0.4134,0.0837,0.5363,0.0275,0.5423,0.2018,0.4728,0.4132,0.1846,0.4502,0.2143,0.4941,0.1944,0.2148,0.6103,0.1697,0.6911,0.753,0.0914,0.8521,-0.019,0.0185,0.0419,0.1373,0.2975");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3_cP60_212_bcd_ace"){
        vparameters.push_back("8.4,0.0,0.871,0.375,0.378,0.129,0.629");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B2_cP20_213_d_c"){
        vparameters.push_back("6.9352,0.07378,0.2051");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC3_cP24_213_c_a_d"){
        vparameters.push_back("6.84,0.061,0.206");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC4_cF24_216_b_a_e"){
        vparameters.push_back("7.5,0.125");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF40_216_de_ce"){
        vparameters.push_back("8.1926,0.1137,0.8804");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB4C8_cF52_216_a_e_2e"){
        vparameters.push_back("9.7294,0.3974,0.6343,0.135");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A13BC18D20E5_cF228_216_dh_b_fh_2eh_ce"){
        vparameters.push_back("13.87,0.825,0.1818,0.1143,0.278,0.0853,0.7667,0.1793,0.5466,0.1385,0.0003");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A16B40C12D6E5_cF316_216_eh_e2g2h_h_f_be"){
        vparameters.push_back("14.886,0.079,0.4262,0.1722,0.1812,0.267,0.515,0.894,-0.085,-0.0729,0.7024,0.6038,0.258,0.8371,0.4913");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A45B11_cF448_216_bd4efg5h_ac2eh"){
        vparameters.push_back("21.699,0.0834,0.9126,0.1636,0.8297,0.4059,0.6618,0.1573,0.0895,0.2958,0.3904,0.4377,0.2627,0.5455,0.6403,0.6728,0.5123,0.9161,0.7637,0.1735,0.0142");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cI16_217_c_c"){
        vparameters.push_back("5.9487,0.3409,0.1657");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B24C13_cI82_217_c_deg_ag"){
        vparameters.push_back("10.1439,0.24271,0.2173,0.0777,0.2123,0.61822,0.14266");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BC4D12E3_cP46_218_d_a_e_i_c"){
        vparameters.push_back("8.882,0.1778,0.139,0.1494,0.4383");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B4C4D4E16F4G3_cP76_218_c_e_e_e_ei_e_d"){
        vparameters.push_back("9.1097,0.2009,0.164,0.2392,0.0995,-0.0331,0.6443,0.6558,-0.0331");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cI24_220_a_b"){
        vparameters.push_back("6.2951");  // (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B12C3_cI76_220_c_e_a"){
        vparameters.push_back("10.2867,0.0849,0.059,0.131,0.2889");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A7B12C19_cI152_220_bc_2d_ace"){
        vparameters.push_back("11.9794,0.0188,0.43519,0.1432,0.1867,0.78672,0.09946,0.30708");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB9C3_cI208_220_c_3e_e"){
        vparameters.push_back("13.63,0.117,0.09,0.11,0.8,0.095,0.141,0.245,0.137,0.096,-0.014,0.34,0.063,0.124");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB12C_cP14_221_a_h_b"){
        vparameters.push_back("3.4708,0.24126");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B40CD12_cP112_224_d_e3k_a_k"){
        vparameters.push_back("11.75,0.3226,0.6463,0.0133,0.3687,0.5462,0.5629,0.2412,0.4636,0.2579");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A5B40CD12_cP116_224_cd_e3k_a_k"){
        vparameters.push_back("12.141,0.83154,0.67011,0.51608,0.87272,0.04157,0.43783,0.74176,-0.0445,0.75577");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A27B52CD12_cP184_224_dl_eh3k_a_k"){
        vparameters.push_back("12.506,0.32273,0.34752,0.65683,-0.00684,0.37155,0.52447,0.55441,0.23276,0.45688,0.25822,0.24525,0.61259,0.68288");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BC6D_cF40_225_c_a_e_b"){
        vparameters.push_back("8.19849,0.2654");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A10B_cF44_225_cf_b"){
        vparameters.push_back("4.78,0.12");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A9B8_cF68_225_af_ce"){
        vparameters.push_back("9.928,0.7591,0.374");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="ABC6D8E2_cF72_225_b_a_e_f_c"){
        vparameters.push_back("10.08,0.226,0.164");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B9CD2E6_cF96_225_e_bf_a_c_e"){
        vparameters.push_back("10.12,0.20257,0.31126,0.16667");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB30C16D3_cF200_225_a_ej_2f_bc"){
        vparameters.push_back("8.94011,0.1969,0.307,0.557,0.0536,0.1923");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB_cF32_227_c_d"){
        vparameters.push_back("7.23");  // (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B_cF96_227_abf_cd"){
        vparameters.push_back("10.24,0.395");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A11B4_cF120_227_acdf_e"){
        vparameters.push_back("8.2376,0.2552,0.368");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A16B2C_cF152_227_eg_d_a"){
        vparameters.push_back("6.71851,0.2915,0.56034,0.8706");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A18B2C3_cF184_227_fg_d_ac"){
        vparameters.push_back("14.53,0.7657,0.5584,0.3252");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A22B_cF184_227_cdfg_a"){
        vparameters.push_back("14.101,0.765,0.562,0.3195");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD3E6_cF208_227_e_c_d_f_g"){
        vparameters.push_back("14.05,0.715,0.35,0.75,0.65");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A4B2C6D16E_cF232_227_e_d_f_eg_a"){
        vparameters.push_back("13.89999,0.735,0.06,0.35,0.765,0.665");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A29B40CD12_cF656_227_ae2fg_e3g_b_g"){
        vparameters.push_back("23.28001,-0.05284,0.41753,-0.01933,0.85395,0.71864,0.85825,0.83333,0.25301,-0.061,0.52706,0.5378,0.3707,-0.01783,0.38144");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A21B_cI44_229_bdh_a"){
        vparameters.push_back("5.034,0.375");  // 001, (part 3)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2B3C12D12_cI232_230_a_c_h_h"){
        vparameters.push_back("12.575,-0.091,0.054,0.163,-0.0282,0.0531,0.1396");  // 001, (part 3)
      }
      // -------------------------------------------------------------------------
      // carbo-nitride prototypes (from DX) //DX20210728
      // -------------------------------------------------------------------------
      // quaternaries
      // ---------------------------------------------------------------------------
      if(anrl_label=="A2BCD4_mP16_11_2e_e_e_4e"){
        vparameters.push_back("10.13,0.37429418,0.53543929,101.35,0.8907,0.574,0.5259,0.7797,0.19343,0.98319,0.6426,0.4055,0.8415,0.3498,0.9342,0.8038,0.4048,0.7576,0.6509,0.7929");  // 001, ternary metal-carbo-nitride prototype (ICSD #420440)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="AB2C2D_mC24_15_a_f_f_e"){
        vparameters.push_back("6.572,0.56451613,2.7888162,108.9486,0.7889,0.6658,0.6423,0.6089,0.761,0.7178,0.6718");  // 001, ternary metal-carbo-nitride prototype (ICSD #65697)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC6D2_hP15_147_g_a_g_d"){
        vparameters.push_back("7.1346,0.76431755,0.60052,0.06546,0.25027,0.19568,0.30145,0.90535,0.30707");  // 001, ternary metal-carbo-nitride prototype (ICSD #51493)
        vparameters.push_back("7.1805,0.74700926,0.57625,0.74903,0.83501,0.19801,0.13274,0.4025,0.31421");  // 002, ternary metal-carbo-nitride prototype (ICSD #51494)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2C6D_hP15_147_g_d_g_a"){
        vparameters.push_back("6.3757,1.033361,0.5068,0.032,0.266,0.164,0.003,0.399,0.259");  // 001, ternary metal-carbo-nitride prototype (ICSD #417821)
        vparameters.push_back("6.3797,1.0313808,0.499,0.997,0.253,0.169,0.988,0.379,0.284");  // 002, ternary metal-carbo-nitride prototype (ICSD #417822)
        vparameters.push_back("6.1341,1.0860273,0.503,0.036,0.289,0.185,0,0.388,0.308");  // 003, ternary metal-carbo-nitride prototype (ICSD #417823)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2CD6_hP15_147_g_d_a_g"){
        vparameters.push_back("6.3205,1.0027529,0.4965,0.016,0.249,0.152,0.03,0.381,0.288");  // 001, ternary metal-carbo-nitride prototype (ICSD #417826)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3B6CD6_hP16_162_g_k_a_k"){
        vparameters.push_back("7.03,1.013798,0.216,0.159,0.346,0.255");  // 001, ternary metal-carbo-nitride prototype (ICSD #16959)
        vparameters.push_back("6.7172,1.1194545,0.22468,0.84235,0.34271,0.73089");  // 002, ternary metal-carbo-nitride prototype (ICSD #173552)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6B2CD6_cF60_216_f_bc_a_f"){
        vparameters.push_back("9.98,0.81,0.69");  // 001, ternary metal-carbo-nitride prototype (ICSD #28653)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A6BC6D_cF56_225_e_a_e_b"){
        vparameters.push_back("10.911,0.3102,0.2081");  // 001, ternary metal-carbo-nitride prototype (ICSD #6093)
      }
      // ---------------------------------------------------------------------------
      if(anrl_label=="A3BCD3_cF64_225_e_c_ab_e"){
        vparameters.push_back("10.6061,0.68214,0.79359");  // 001, ternary metal-carbo-nitride prototype (ICSD #248046)
      }
#endif
    }
    if(library=="" && !keep_original_lattice_parameter){
      //DX 20190314 - loop over parameters - START
      for(uint p=0;p<vparameters.size();p++){
        aurostd::string2tokens(vparameters[p],tokens,",");
        tokens[0]="-1";
        vparameters[p]=aurostd::joinWDelimiter(tokens,",");
      }
      //DX 20190314 - loop over parameters - END
    }
    if(library != "all"){
      if(library=="" && choice!=-1){
        if((uint)choice<vparameters.size()){
          vector<string> tmp; tmp.push_back(vparameters[choice]);
          vparameters.clear();
          vparameters=tmp;
        }
        else{
          message << "anrl::getANRLParameters(): ERROR - " << anrl_label << " does not have more than " << vparameters.size() << " choice(s):";
          for(uint i=0;i<vparameters.size();i++){
            message << "  " << anrl_label << "-" << std::setw(3) << std::setfill('0') << i+1 << " : " << vparameters[i] << endl;
          }
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _VALUE_RANGE_); //DX 20191118 - exit to throw
        }
      }
      else if((library=="" && choice==-1) || (vparameters.size() && (library=="part1" || library=="part2" || library=="misc"))){
        // -------------------------------------------------------------------------
        // we do not automatically assign the enumeration suffix to the label
        // (even if there there is currently only one prototype with a label)
        // otherwise, the command without the suffix may break in the future
        // (cannot auto-assign when more than one choice, e.g., -001 or -002)
        message << "anrl::getANRLParameters(): ERROR - " << anrl_label << " has " << vparameters.size() << " preset parameter set(s): " << endl;
        for(uint i=0;i<vparameters.size();i++){
          message << "  " << anrl_label << "-" << std::setw(3) << std::setfill('0') << i+1 << " : " << vparameters[i] << endl;
        }
        message << "Rerun command and specify the parameters or the preset suffix, e.g., aflow --proto=" << anrl_label << "-" << std::setw(3) << std::setfill('0') << 1; //DX 20190826 - changed "./aflow" to "aflow"
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _VALUE_ERROR_); //DX 20191118 - exit to throw
      }
    }
    return vparameters;  
  }
}

#endif // AFLOW_REMOVE_GREP // _AFLOW_ANRL_LIST_CPP // AFLOW_REMOVE_GREP

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************

