/* Produced by CVXGEN, 2021-08-10 14:54:42 -0400.  */
/* CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2017 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Filename: matrix_support.c. */
/* Description: Support functions for matrix multiplication and vector filling. */
#include "solver.h"
void multbymA(double *lhs, double *rhs) {
}
void multbymAT(double *lhs, double *rhs) {
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
}
void multbymG(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(params.A[0])-rhs[1]*(params.A[60])-rhs[2]*(params.A[120])-rhs[3]*(params.A[180])-rhs[4]*(params.A[240])-rhs[5]*(params.A[300])-rhs[6]*(params.A[360])-rhs[7]*(params.A[420])-rhs[8]*(params.A[480])-rhs[9]*(params.A[540]);
  lhs[1] = -rhs[0]*(params.A[1])-rhs[1]*(params.A[61])-rhs[2]*(params.A[121])-rhs[3]*(params.A[181])-rhs[4]*(params.A[241])-rhs[5]*(params.A[301])-rhs[6]*(params.A[361])-rhs[7]*(params.A[421])-rhs[8]*(params.A[481])-rhs[9]*(params.A[541]);
  lhs[2] = -rhs[0]*(params.A[2])-rhs[1]*(params.A[62])-rhs[2]*(params.A[122])-rhs[3]*(params.A[182])-rhs[4]*(params.A[242])-rhs[5]*(params.A[302])-rhs[6]*(params.A[362])-rhs[7]*(params.A[422])-rhs[8]*(params.A[482])-rhs[9]*(params.A[542]);
  lhs[3] = -rhs[0]*(params.A[3])-rhs[1]*(params.A[63])-rhs[2]*(params.A[123])-rhs[3]*(params.A[183])-rhs[4]*(params.A[243])-rhs[5]*(params.A[303])-rhs[6]*(params.A[363])-rhs[7]*(params.A[423])-rhs[8]*(params.A[483])-rhs[9]*(params.A[543]);
  lhs[4] = -rhs[0]*(params.A[4])-rhs[1]*(params.A[64])-rhs[2]*(params.A[124])-rhs[3]*(params.A[184])-rhs[4]*(params.A[244])-rhs[5]*(params.A[304])-rhs[6]*(params.A[364])-rhs[7]*(params.A[424])-rhs[8]*(params.A[484])-rhs[9]*(params.A[544]);
  lhs[5] = -rhs[0]*(params.A[5])-rhs[1]*(params.A[65])-rhs[2]*(params.A[125])-rhs[3]*(params.A[185])-rhs[4]*(params.A[245])-rhs[5]*(params.A[305])-rhs[6]*(params.A[365])-rhs[7]*(params.A[425])-rhs[8]*(params.A[485])-rhs[9]*(params.A[545]);
  lhs[6] = -rhs[0]*(params.A[6])-rhs[1]*(params.A[66])-rhs[2]*(params.A[126])-rhs[3]*(params.A[186])-rhs[4]*(params.A[246])-rhs[5]*(params.A[306])-rhs[6]*(params.A[366])-rhs[7]*(params.A[426])-rhs[8]*(params.A[486])-rhs[9]*(params.A[546]);
  lhs[7] = -rhs[0]*(params.A[7])-rhs[1]*(params.A[67])-rhs[2]*(params.A[127])-rhs[3]*(params.A[187])-rhs[4]*(params.A[247])-rhs[5]*(params.A[307])-rhs[6]*(params.A[367])-rhs[7]*(params.A[427])-rhs[8]*(params.A[487])-rhs[9]*(params.A[547]);
  lhs[8] = -rhs[0]*(params.A[8])-rhs[1]*(params.A[68])-rhs[2]*(params.A[128])-rhs[3]*(params.A[188])-rhs[4]*(params.A[248])-rhs[5]*(params.A[308])-rhs[6]*(params.A[368])-rhs[7]*(params.A[428])-rhs[8]*(params.A[488])-rhs[9]*(params.A[548]);
  lhs[9] = -rhs[0]*(params.A[9])-rhs[1]*(params.A[69])-rhs[2]*(params.A[129])-rhs[3]*(params.A[189])-rhs[4]*(params.A[249])-rhs[5]*(params.A[309])-rhs[6]*(params.A[369])-rhs[7]*(params.A[429])-rhs[8]*(params.A[489])-rhs[9]*(params.A[549]);
  lhs[10] = -rhs[0]*(params.A[10])-rhs[1]*(params.A[70])-rhs[2]*(params.A[130])-rhs[3]*(params.A[190])-rhs[4]*(params.A[250])-rhs[5]*(params.A[310])-rhs[6]*(params.A[370])-rhs[7]*(params.A[430])-rhs[8]*(params.A[490])-rhs[9]*(params.A[550]);
  lhs[11] = -rhs[0]*(params.A[11])-rhs[1]*(params.A[71])-rhs[2]*(params.A[131])-rhs[3]*(params.A[191])-rhs[4]*(params.A[251])-rhs[5]*(params.A[311])-rhs[6]*(params.A[371])-rhs[7]*(params.A[431])-rhs[8]*(params.A[491])-rhs[9]*(params.A[551]);
  lhs[12] = -rhs[0]*(params.A[12])-rhs[1]*(params.A[72])-rhs[2]*(params.A[132])-rhs[3]*(params.A[192])-rhs[4]*(params.A[252])-rhs[5]*(params.A[312])-rhs[6]*(params.A[372])-rhs[7]*(params.A[432])-rhs[8]*(params.A[492])-rhs[9]*(params.A[552]);
  lhs[13] = -rhs[0]*(params.A[13])-rhs[1]*(params.A[73])-rhs[2]*(params.A[133])-rhs[3]*(params.A[193])-rhs[4]*(params.A[253])-rhs[5]*(params.A[313])-rhs[6]*(params.A[373])-rhs[7]*(params.A[433])-rhs[8]*(params.A[493])-rhs[9]*(params.A[553]);
  lhs[14] = -rhs[0]*(params.A[14])-rhs[1]*(params.A[74])-rhs[2]*(params.A[134])-rhs[3]*(params.A[194])-rhs[4]*(params.A[254])-rhs[5]*(params.A[314])-rhs[6]*(params.A[374])-rhs[7]*(params.A[434])-rhs[8]*(params.A[494])-rhs[9]*(params.A[554]);
  lhs[15] = -rhs[0]*(params.A[15])-rhs[1]*(params.A[75])-rhs[2]*(params.A[135])-rhs[3]*(params.A[195])-rhs[4]*(params.A[255])-rhs[5]*(params.A[315])-rhs[6]*(params.A[375])-rhs[7]*(params.A[435])-rhs[8]*(params.A[495])-rhs[9]*(params.A[555]);
  lhs[16] = -rhs[0]*(params.A[16])-rhs[1]*(params.A[76])-rhs[2]*(params.A[136])-rhs[3]*(params.A[196])-rhs[4]*(params.A[256])-rhs[5]*(params.A[316])-rhs[6]*(params.A[376])-rhs[7]*(params.A[436])-rhs[8]*(params.A[496])-rhs[9]*(params.A[556]);
  lhs[17] = -rhs[0]*(params.A[17])-rhs[1]*(params.A[77])-rhs[2]*(params.A[137])-rhs[3]*(params.A[197])-rhs[4]*(params.A[257])-rhs[5]*(params.A[317])-rhs[6]*(params.A[377])-rhs[7]*(params.A[437])-rhs[8]*(params.A[497])-rhs[9]*(params.A[557]);
  lhs[18] = -rhs[0]*(params.A[18])-rhs[1]*(params.A[78])-rhs[2]*(params.A[138])-rhs[3]*(params.A[198])-rhs[4]*(params.A[258])-rhs[5]*(params.A[318])-rhs[6]*(params.A[378])-rhs[7]*(params.A[438])-rhs[8]*(params.A[498])-rhs[9]*(params.A[558]);
  lhs[19] = -rhs[0]*(params.A[19])-rhs[1]*(params.A[79])-rhs[2]*(params.A[139])-rhs[3]*(params.A[199])-rhs[4]*(params.A[259])-rhs[5]*(params.A[319])-rhs[6]*(params.A[379])-rhs[7]*(params.A[439])-rhs[8]*(params.A[499])-rhs[9]*(params.A[559]);
  lhs[20] = -rhs[0]*(params.A[20])-rhs[1]*(params.A[80])-rhs[2]*(params.A[140])-rhs[3]*(params.A[200])-rhs[4]*(params.A[260])-rhs[5]*(params.A[320])-rhs[6]*(params.A[380])-rhs[7]*(params.A[440])-rhs[8]*(params.A[500])-rhs[9]*(params.A[560]);
  lhs[21] = -rhs[0]*(params.A[21])-rhs[1]*(params.A[81])-rhs[2]*(params.A[141])-rhs[3]*(params.A[201])-rhs[4]*(params.A[261])-rhs[5]*(params.A[321])-rhs[6]*(params.A[381])-rhs[7]*(params.A[441])-rhs[8]*(params.A[501])-rhs[9]*(params.A[561]);
  lhs[22] = -rhs[0]*(params.A[22])-rhs[1]*(params.A[82])-rhs[2]*(params.A[142])-rhs[3]*(params.A[202])-rhs[4]*(params.A[262])-rhs[5]*(params.A[322])-rhs[6]*(params.A[382])-rhs[7]*(params.A[442])-rhs[8]*(params.A[502])-rhs[9]*(params.A[562]);
  lhs[23] = -rhs[0]*(params.A[23])-rhs[1]*(params.A[83])-rhs[2]*(params.A[143])-rhs[3]*(params.A[203])-rhs[4]*(params.A[263])-rhs[5]*(params.A[323])-rhs[6]*(params.A[383])-rhs[7]*(params.A[443])-rhs[8]*(params.A[503])-rhs[9]*(params.A[563]);
  lhs[24] = -rhs[0]*(params.A[24])-rhs[1]*(params.A[84])-rhs[2]*(params.A[144])-rhs[3]*(params.A[204])-rhs[4]*(params.A[264])-rhs[5]*(params.A[324])-rhs[6]*(params.A[384])-rhs[7]*(params.A[444])-rhs[8]*(params.A[504])-rhs[9]*(params.A[564]);
  lhs[25] = -rhs[0]*(params.A[25])-rhs[1]*(params.A[85])-rhs[2]*(params.A[145])-rhs[3]*(params.A[205])-rhs[4]*(params.A[265])-rhs[5]*(params.A[325])-rhs[6]*(params.A[385])-rhs[7]*(params.A[445])-rhs[8]*(params.A[505])-rhs[9]*(params.A[565]);
  lhs[26] = -rhs[0]*(params.A[26])-rhs[1]*(params.A[86])-rhs[2]*(params.A[146])-rhs[3]*(params.A[206])-rhs[4]*(params.A[266])-rhs[5]*(params.A[326])-rhs[6]*(params.A[386])-rhs[7]*(params.A[446])-rhs[8]*(params.A[506])-rhs[9]*(params.A[566]);
  lhs[27] = -rhs[0]*(params.A[27])-rhs[1]*(params.A[87])-rhs[2]*(params.A[147])-rhs[3]*(params.A[207])-rhs[4]*(params.A[267])-rhs[5]*(params.A[327])-rhs[6]*(params.A[387])-rhs[7]*(params.A[447])-rhs[8]*(params.A[507])-rhs[9]*(params.A[567]);
  lhs[28] = -rhs[0]*(params.A[28])-rhs[1]*(params.A[88])-rhs[2]*(params.A[148])-rhs[3]*(params.A[208])-rhs[4]*(params.A[268])-rhs[5]*(params.A[328])-rhs[6]*(params.A[388])-rhs[7]*(params.A[448])-rhs[8]*(params.A[508])-rhs[9]*(params.A[568]);
  lhs[29] = -rhs[0]*(params.A[29])-rhs[1]*(params.A[89])-rhs[2]*(params.A[149])-rhs[3]*(params.A[209])-rhs[4]*(params.A[269])-rhs[5]*(params.A[329])-rhs[6]*(params.A[389])-rhs[7]*(params.A[449])-rhs[8]*(params.A[509])-rhs[9]*(params.A[569]);
  lhs[30] = -rhs[0]*(params.A[30])-rhs[1]*(params.A[90])-rhs[2]*(params.A[150])-rhs[3]*(params.A[210])-rhs[4]*(params.A[270])-rhs[5]*(params.A[330])-rhs[6]*(params.A[390])-rhs[7]*(params.A[450])-rhs[8]*(params.A[510])-rhs[9]*(params.A[570]);
  lhs[31] = -rhs[0]*(params.A[31])-rhs[1]*(params.A[91])-rhs[2]*(params.A[151])-rhs[3]*(params.A[211])-rhs[4]*(params.A[271])-rhs[5]*(params.A[331])-rhs[6]*(params.A[391])-rhs[7]*(params.A[451])-rhs[8]*(params.A[511])-rhs[9]*(params.A[571]);
  lhs[32] = -rhs[0]*(params.A[32])-rhs[1]*(params.A[92])-rhs[2]*(params.A[152])-rhs[3]*(params.A[212])-rhs[4]*(params.A[272])-rhs[5]*(params.A[332])-rhs[6]*(params.A[392])-rhs[7]*(params.A[452])-rhs[8]*(params.A[512])-rhs[9]*(params.A[572]);
  lhs[33] = -rhs[0]*(params.A[33])-rhs[1]*(params.A[93])-rhs[2]*(params.A[153])-rhs[3]*(params.A[213])-rhs[4]*(params.A[273])-rhs[5]*(params.A[333])-rhs[6]*(params.A[393])-rhs[7]*(params.A[453])-rhs[8]*(params.A[513])-rhs[9]*(params.A[573]);
  lhs[34] = -rhs[0]*(params.A[34])-rhs[1]*(params.A[94])-rhs[2]*(params.A[154])-rhs[3]*(params.A[214])-rhs[4]*(params.A[274])-rhs[5]*(params.A[334])-rhs[6]*(params.A[394])-rhs[7]*(params.A[454])-rhs[8]*(params.A[514])-rhs[9]*(params.A[574]);
  lhs[35] = -rhs[0]*(params.A[35])-rhs[1]*(params.A[95])-rhs[2]*(params.A[155])-rhs[3]*(params.A[215])-rhs[4]*(params.A[275])-rhs[5]*(params.A[335])-rhs[6]*(params.A[395])-rhs[7]*(params.A[455])-rhs[8]*(params.A[515])-rhs[9]*(params.A[575]);
  lhs[36] = -rhs[0]*(params.A[36])-rhs[1]*(params.A[96])-rhs[2]*(params.A[156])-rhs[3]*(params.A[216])-rhs[4]*(params.A[276])-rhs[5]*(params.A[336])-rhs[6]*(params.A[396])-rhs[7]*(params.A[456])-rhs[8]*(params.A[516])-rhs[9]*(params.A[576]);
  lhs[37] = -rhs[0]*(params.A[37])-rhs[1]*(params.A[97])-rhs[2]*(params.A[157])-rhs[3]*(params.A[217])-rhs[4]*(params.A[277])-rhs[5]*(params.A[337])-rhs[6]*(params.A[397])-rhs[7]*(params.A[457])-rhs[8]*(params.A[517])-rhs[9]*(params.A[577]);
  lhs[38] = -rhs[0]*(params.A[38])-rhs[1]*(params.A[98])-rhs[2]*(params.A[158])-rhs[3]*(params.A[218])-rhs[4]*(params.A[278])-rhs[5]*(params.A[338])-rhs[6]*(params.A[398])-rhs[7]*(params.A[458])-rhs[8]*(params.A[518])-rhs[9]*(params.A[578]);
  lhs[39] = -rhs[0]*(params.A[39])-rhs[1]*(params.A[99])-rhs[2]*(params.A[159])-rhs[3]*(params.A[219])-rhs[4]*(params.A[279])-rhs[5]*(params.A[339])-rhs[6]*(params.A[399])-rhs[7]*(params.A[459])-rhs[8]*(params.A[519])-rhs[9]*(params.A[579]);
  lhs[40] = -rhs[0]*(params.A[40])-rhs[1]*(params.A[100])-rhs[2]*(params.A[160])-rhs[3]*(params.A[220])-rhs[4]*(params.A[280])-rhs[5]*(params.A[340])-rhs[6]*(params.A[400])-rhs[7]*(params.A[460])-rhs[8]*(params.A[520])-rhs[9]*(params.A[580]);
  lhs[41] = -rhs[0]*(params.A[41])-rhs[1]*(params.A[101])-rhs[2]*(params.A[161])-rhs[3]*(params.A[221])-rhs[4]*(params.A[281])-rhs[5]*(params.A[341])-rhs[6]*(params.A[401])-rhs[7]*(params.A[461])-rhs[8]*(params.A[521])-rhs[9]*(params.A[581]);
  lhs[42] = -rhs[0]*(params.A[42])-rhs[1]*(params.A[102])-rhs[2]*(params.A[162])-rhs[3]*(params.A[222])-rhs[4]*(params.A[282])-rhs[5]*(params.A[342])-rhs[6]*(params.A[402])-rhs[7]*(params.A[462])-rhs[8]*(params.A[522])-rhs[9]*(params.A[582]);
  lhs[43] = -rhs[0]*(params.A[43])-rhs[1]*(params.A[103])-rhs[2]*(params.A[163])-rhs[3]*(params.A[223])-rhs[4]*(params.A[283])-rhs[5]*(params.A[343])-rhs[6]*(params.A[403])-rhs[7]*(params.A[463])-rhs[8]*(params.A[523])-rhs[9]*(params.A[583]);
  lhs[44] = -rhs[0]*(params.A[44])-rhs[1]*(params.A[104])-rhs[2]*(params.A[164])-rhs[3]*(params.A[224])-rhs[4]*(params.A[284])-rhs[5]*(params.A[344])-rhs[6]*(params.A[404])-rhs[7]*(params.A[464])-rhs[8]*(params.A[524])-rhs[9]*(params.A[584]);
  lhs[45] = -rhs[0]*(params.A[45])-rhs[1]*(params.A[105])-rhs[2]*(params.A[165])-rhs[3]*(params.A[225])-rhs[4]*(params.A[285])-rhs[5]*(params.A[345])-rhs[6]*(params.A[405])-rhs[7]*(params.A[465])-rhs[8]*(params.A[525])-rhs[9]*(params.A[585]);
  lhs[46] = -rhs[0]*(params.A[46])-rhs[1]*(params.A[106])-rhs[2]*(params.A[166])-rhs[3]*(params.A[226])-rhs[4]*(params.A[286])-rhs[5]*(params.A[346])-rhs[6]*(params.A[406])-rhs[7]*(params.A[466])-rhs[8]*(params.A[526])-rhs[9]*(params.A[586]);
  lhs[47] = -rhs[0]*(params.A[47])-rhs[1]*(params.A[107])-rhs[2]*(params.A[167])-rhs[3]*(params.A[227])-rhs[4]*(params.A[287])-rhs[5]*(params.A[347])-rhs[6]*(params.A[407])-rhs[7]*(params.A[467])-rhs[8]*(params.A[527])-rhs[9]*(params.A[587]);
  lhs[48] = -rhs[0]*(params.A[48])-rhs[1]*(params.A[108])-rhs[2]*(params.A[168])-rhs[3]*(params.A[228])-rhs[4]*(params.A[288])-rhs[5]*(params.A[348])-rhs[6]*(params.A[408])-rhs[7]*(params.A[468])-rhs[8]*(params.A[528])-rhs[9]*(params.A[588]);
  lhs[49] = -rhs[0]*(params.A[49])-rhs[1]*(params.A[109])-rhs[2]*(params.A[169])-rhs[3]*(params.A[229])-rhs[4]*(params.A[289])-rhs[5]*(params.A[349])-rhs[6]*(params.A[409])-rhs[7]*(params.A[469])-rhs[8]*(params.A[529])-rhs[9]*(params.A[589]);
  lhs[50] = -rhs[0]*(params.A[50])-rhs[1]*(params.A[110])-rhs[2]*(params.A[170])-rhs[3]*(params.A[230])-rhs[4]*(params.A[290])-rhs[5]*(params.A[350])-rhs[6]*(params.A[410])-rhs[7]*(params.A[470])-rhs[8]*(params.A[530])-rhs[9]*(params.A[590]);
  lhs[51] = -rhs[0]*(params.A[51])-rhs[1]*(params.A[111])-rhs[2]*(params.A[171])-rhs[3]*(params.A[231])-rhs[4]*(params.A[291])-rhs[5]*(params.A[351])-rhs[6]*(params.A[411])-rhs[7]*(params.A[471])-rhs[8]*(params.A[531])-rhs[9]*(params.A[591]);
  lhs[52] = -rhs[0]*(params.A[52])-rhs[1]*(params.A[112])-rhs[2]*(params.A[172])-rhs[3]*(params.A[232])-rhs[4]*(params.A[292])-rhs[5]*(params.A[352])-rhs[6]*(params.A[412])-rhs[7]*(params.A[472])-rhs[8]*(params.A[532])-rhs[9]*(params.A[592]);
  lhs[53] = -rhs[0]*(params.A[53])-rhs[1]*(params.A[113])-rhs[2]*(params.A[173])-rhs[3]*(params.A[233])-rhs[4]*(params.A[293])-rhs[5]*(params.A[353])-rhs[6]*(params.A[413])-rhs[7]*(params.A[473])-rhs[8]*(params.A[533])-rhs[9]*(params.A[593]);
  lhs[54] = -rhs[0]*(params.A[54])-rhs[1]*(params.A[114])-rhs[2]*(params.A[174])-rhs[3]*(params.A[234])-rhs[4]*(params.A[294])-rhs[5]*(params.A[354])-rhs[6]*(params.A[414])-rhs[7]*(params.A[474])-rhs[8]*(params.A[534])-rhs[9]*(params.A[594]);
  lhs[55] = -rhs[0]*(params.A[55])-rhs[1]*(params.A[115])-rhs[2]*(params.A[175])-rhs[3]*(params.A[235])-rhs[4]*(params.A[295])-rhs[5]*(params.A[355])-rhs[6]*(params.A[415])-rhs[7]*(params.A[475])-rhs[8]*(params.A[535])-rhs[9]*(params.A[595]);
  lhs[56] = -rhs[0]*(params.A[56])-rhs[1]*(params.A[116])-rhs[2]*(params.A[176])-rhs[3]*(params.A[236])-rhs[4]*(params.A[296])-rhs[5]*(params.A[356])-rhs[6]*(params.A[416])-rhs[7]*(params.A[476])-rhs[8]*(params.A[536])-rhs[9]*(params.A[596]);
  lhs[57] = -rhs[0]*(params.A[57])-rhs[1]*(params.A[117])-rhs[2]*(params.A[177])-rhs[3]*(params.A[237])-rhs[4]*(params.A[297])-rhs[5]*(params.A[357])-rhs[6]*(params.A[417])-rhs[7]*(params.A[477])-rhs[8]*(params.A[537])-rhs[9]*(params.A[597]);
  lhs[58] = -rhs[0]*(params.A[58])-rhs[1]*(params.A[118])-rhs[2]*(params.A[178])-rhs[3]*(params.A[238])-rhs[4]*(params.A[298])-rhs[5]*(params.A[358])-rhs[6]*(params.A[418])-rhs[7]*(params.A[478])-rhs[8]*(params.A[538])-rhs[9]*(params.A[598]);
  lhs[59] = -rhs[0]*(params.A[59])-rhs[1]*(params.A[119])-rhs[2]*(params.A[179])-rhs[3]*(params.A[239])-rhs[4]*(params.A[299])-rhs[5]*(params.A[359])-rhs[6]*(params.A[419])-rhs[7]*(params.A[479])-rhs[8]*(params.A[539])-rhs[9]*(params.A[599]);
}
void multbymGT(double *lhs, double *rhs) {
  lhs[0] = -rhs[0]*(params.A[0])-rhs[1]*(params.A[1])-rhs[2]*(params.A[2])-rhs[3]*(params.A[3])-rhs[4]*(params.A[4])-rhs[5]*(params.A[5])-rhs[6]*(params.A[6])-rhs[7]*(params.A[7])-rhs[8]*(params.A[8])-rhs[9]*(params.A[9])-rhs[10]*(params.A[10])-rhs[11]*(params.A[11])-rhs[12]*(params.A[12])-rhs[13]*(params.A[13])-rhs[14]*(params.A[14])-rhs[15]*(params.A[15])-rhs[16]*(params.A[16])-rhs[17]*(params.A[17])-rhs[18]*(params.A[18])-rhs[19]*(params.A[19])-rhs[20]*(params.A[20])-rhs[21]*(params.A[21])-rhs[22]*(params.A[22])-rhs[23]*(params.A[23])-rhs[24]*(params.A[24])-rhs[25]*(params.A[25])-rhs[26]*(params.A[26])-rhs[27]*(params.A[27])-rhs[28]*(params.A[28])-rhs[29]*(params.A[29])-rhs[30]*(params.A[30])-rhs[31]*(params.A[31])-rhs[32]*(params.A[32])-rhs[33]*(params.A[33])-rhs[34]*(params.A[34])-rhs[35]*(params.A[35])-rhs[36]*(params.A[36])-rhs[37]*(params.A[37])-rhs[38]*(params.A[38])-rhs[39]*(params.A[39])-rhs[40]*(params.A[40])-rhs[41]*(params.A[41])-rhs[42]*(params.A[42])-rhs[43]*(params.A[43])-rhs[44]*(params.A[44])-rhs[45]*(params.A[45])-rhs[46]*(params.A[46])-rhs[47]*(params.A[47])-rhs[48]*(params.A[48])-rhs[49]*(params.A[49])-rhs[50]*(params.A[50])-rhs[51]*(params.A[51])-rhs[52]*(params.A[52])-rhs[53]*(params.A[53])-rhs[54]*(params.A[54])-rhs[55]*(params.A[55])-rhs[56]*(params.A[56])-rhs[57]*(params.A[57])-rhs[58]*(params.A[58])-rhs[59]*(params.A[59]);
  lhs[1] = -rhs[0]*(params.A[60])-rhs[1]*(params.A[61])-rhs[2]*(params.A[62])-rhs[3]*(params.A[63])-rhs[4]*(params.A[64])-rhs[5]*(params.A[65])-rhs[6]*(params.A[66])-rhs[7]*(params.A[67])-rhs[8]*(params.A[68])-rhs[9]*(params.A[69])-rhs[10]*(params.A[70])-rhs[11]*(params.A[71])-rhs[12]*(params.A[72])-rhs[13]*(params.A[73])-rhs[14]*(params.A[74])-rhs[15]*(params.A[75])-rhs[16]*(params.A[76])-rhs[17]*(params.A[77])-rhs[18]*(params.A[78])-rhs[19]*(params.A[79])-rhs[20]*(params.A[80])-rhs[21]*(params.A[81])-rhs[22]*(params.A[82])-rhs[23]*(params.A[83])-rhs[24]*(params.A[84])-rhs[25]*(params.A[85])-rhs[26]*(params.A[86])-rhs[27]*(params.A[87])-rhs[28]*(params.A[88])-rhs[29]*(params.A[89])-rhs[30]*(params.A[90])-rhs[31]*(params.A[91])-rhs[32]*(params.A[92])-rhs[33]*(params.A[93])-rhs[34]*(params.A[94])-rhs[35]*(params.A[95])-rhs[36]*(params.A[96])-rhs[37]*(params.A[97])-rhs[38]*(params.A[98])-rhs[39]*(params.A[99])-rhs[40]*(params.A[100])-rhs[41]*(params.A[101])-rhs[42]*(params.A[102])-rhs[43]*(params.A[103])-rhs[44]*(params.A[104])-rhs[45]*(params.A[105])-rhs[46]*(params.A[106])-rhs[47]*(params.A[107])-rhs[48]*(params.A[108])-rhs[49]*(params.A[109])-rhs[50]*(params.A[110])-rhs[51]*(params.A[111])-rhs[52]*(params.A[112])-rhs[53]*(params.A[113])-rhs[54]*(params.A[114])-rhs[55]*(params.A[115])-rhs[56]*(params.A[116])-rhs[57]*(params.A[117])-rhs[58]*(params.A[118])-rhs[59]*(params.A[119]);
  lhs[2] = -rhs[0]*(params.A[120])-rhs[1]*(params.A[121])-rhs[2]*(params.A[122])-rhs[3]*(params.A[123])-rhs[4]*(params.A[124])-rhs[5]*(params.A[125])-rhs[6]*(params.A[126])-rhs[7]*(params.A[127])-rhs[8]*(params.A[128])-rhs[9]*(params.A[129])-rhs[10]*(params.A[130])-rhs[11]*(params.A[131])-rhs[12]*(params.A[132])-rhs[13]*(params.A[133])-rhs[14]*(params.A[134])-rhs[15]*(params.A[135])-rhs[16]*(params.A[136])-rhs[17]*(params.A[137])-rhs[18]*(params.A[138])-rhs[19]*(params.A[139])-rhs[20]*(params.A[140])-rhs[21]*(params.A[141])-rhs[22]*(params.A[142])-rhs[23]*(params.A[143])-rhs[24]*(params.A[144])-rhs[25]*(params.A[145])-rhs[26]*(params.A[146])-rhs[27]*(params.A[147])-rhs[28]*(params.A[148])-rhs[29]*(params.A[149])-rhs[30]*(params.A[150])-rhs[31]*(params.A[151])-rhs[32]*(params.A[152])-rhs[33]*(params.A[153])-rhs[34]*(params.A[154])-rhs[35]*(params.A[155])-rhs[36]*(params.A[156])-rhs[37]*(params.A[157])-rhs[38]*(params.A[158])-rhs[39]*(params.A[159])-rhs[40]*(params.A[160])-rhs[41]*(params.A[161])-rhs[42]*(params.A[162])-rhs[43]*(params.A[163])-rhs[44]*(params.A[164])-rhs[45]*(params.A[165])-rhs[46]*(params.A[166])-rhs[47]*(params.A[167])-rhs[48]*(params.A[168])-rhs[49]*(params.A[169])-rhs[50]*(params.A[170])-rhs[51]*(params.A[171])-rhs[52]*(params.A[172])-rhs[53]*(params.A[173])-rhs[54]*(params.A[174])-rhs[55]*(params.A[175])-rhs[56]*(params.A[176])-rhs[57]*(params.A[177])-rhs[58]*(params.A[178])-rhs[59]*(params.A[179]);
  lhs[3] = -rhs[0]*(params.A[180])-rhs[1]*(params.A[181])-rhs[2]*(params.A[182])-rhs[3]*(params.A[183])-rhs[4]*(params.A[184])-rhs[5]*(params.A[185])-rhs[6]*(params.A[186])-rhs[7]*(params.A[187])-rhs[8]*(params.A[188])-rhs[9]*(params.A[189])-rhs[10]*(params.A[190])-rhs[11]*(params.A[191])-rhs[12]*(params.A[192])-rhs[13]*(params.A[193])-rhs[14]*(params.A[194])-rhs[15]*(params.A[195])-rhs[16]*(params.A[196])-rhs[17]*(params.A[197])-rhs[18]*(params.A[198])-rhs[19]*(params.A[199])-rhs[20]*(params.A[200])-rhs[21]*(params.A[201])-rhs[22]*(params.A[202])-rhs[23]*(params.A[203])-rhs[24]*(params.A[204])-rhs[25]*(params.A[205])-rhs[26]*(params.A[206])-rhs[27]*(params.A[207])-rhs[28]*(params.A[208])-rhs[29]*(params.A[209])-rhs[30]*(params.A[210])-rhs[31]*(params.A[211])-rhs[32]*(params.A[212])-rhs[33]*(params.A[213])-rhs[34]*(params.A[214])-rhs[35]*(params.A[215])-rhs[36]*(params.A[216])-rhs[37]*(params.A[217])-rhs[38]*(params.A[218])-rhs[39]*(params.A[219])-rhs[40]*(params.A[220])-rhs[41]*(params.A[221])-rhs[42]*(params.A[222])-rhs[43]*(params.A[223])-rhs[44]*(params.A[224])-rhs[45]*(params.A[225])-rhs[46]*(params.A[226])-rhs[47]*(params.A[227])-rhs[48]*(params.A[228])-rhs[49]*(params.A[229])-rhs[50]*(params.A[230])-rhs[51]*(params.A[231])-rhs[52]*(params.A[232])-rhs[53]*(params.A[233])-rhs[54]*(params.A[234])-rhs[55]*(params.A[235])-rhs[56]*(params.A[236])-rhs[57]*(params.A[237])-rhs[58]*(params.A[238])-rhs[59]*(params.A[239]);
  lhs[4] = -rhs[0]*(params.A[240])-rhs[1]*(params.A[241])-rhs[2]*(params.A[242])-rhs[3]*(params.A[243])-rhs[4]*(params.A[244])-rhs[5]*(params.A[245])-rhs[6]*(params.A[246])-rhs[7]*(params.A[247])-rhs[8]*(params.A[248])-rhs[9]*(params.A[249])-rhs[10]*(params.A[250])-rhs[11]*(params.A[251])-rhs[12]*(params.A[252])-rhs[13]*(params.A[253])-rhs[14]*(params.A[254])-rhs[15]*(params.A[255])-rhs[16]*(params.A[256])-rhs[17]*(params.A[257])-rhs[18]*(params.A[258])-rhs[19]*(params.A[259])-rhs[20]*(params.A[260])-rhs[21]*(params.A[261])-rhs[22]*(params.A[262])-rhs[23]*(params.A[263])-rhs[24]*(params.A[264])-rhs[25]*(params.A[265])-rhs[26]*(params.A[266])-rhs[27]*(params.A[267])-rhs[28]*(params.A[268])-rhs[29]*(params.A[269])-rhs[30]*(params.A[270])-rhs[31]*(params.A[271])-rhs[32]*(params.A[272])-rhs[33]*(params.A[273])-rhs[34]*(params.A[274])-rhs[35]*(params.A[275])-rhs[36]*(params.A[276])-rhs[37]*(params.A[277])-rhs[38]*(params.A[278])-rhs[39]*(params.A[279])-rhs[40]*(params.A[280])-rhs[41]*(params.A[281])-rhs[42]*(params.A[282])-rhs[43]*(params.A[283])-rhs[44]*(params.A[284])-rhs[45]*(params.A[285])-rhs[46]*(params.A[286])-rhs[47]*(params.A[287])-rhs[48]*(params.A[288])-rhs[49]*(params.A[289])-rhs[50]*(params.A[290])-rhs[51]*(params.A[291])-rhs[52]*(params.A[292])-rhs[53]*(params.A[293])-rhs[54]*(params.A[294])-rhs[55]*(params.A[295])-rhs[56]*(params.A[296])-rhs[57]*(params.A[297])-rhs[58]*(params.A[298])-rhs[59]*(params.A[299]);
  lhs[5] = -rhs[0]*(params.A[300])-rhs[1]*(params.A[301])-rhs[2]*(params.A[302])-rhs[3]*(params.A[303])-rhs[4]*(params.A[304])-rhs[5]*(params.A[305])-rhs[6]*(params.A[306])-rhs[7]*(params.A[307])-rhs[8]*(params.A[308])-rhs[9]*(params.A[309])-rhs[10]*(params.A[310])-rhs[11]*(params.A[311])-rhs[12]*(params.A[312])-rhs[13]*(params.A[313])-rhs[14]*(params.A[314])-rhs[15]*(params.A[315])-rhs[16]*(params.A[316])-rhs[17]*(params.A[317])-rhs[18]*(params.A[318])-rhs[19]*(params.A[319])-rhs[20]*(params.A[320])-rhs[21]*(params.A[321])-rhs[22]*(params.A[322])-rhs[23]*(params.A[323])-rhs[24]*(params.A[324])-rhs[25]*(params.A[325])-rhs[26]*(params.A[326])-rhs[27]*(params.A[327])-rhs[28]*(params.A[328])-rhs[29]*(params.A[329])-rhs[30]*(params.A[330])-rhs[31]*(params.A[331])-rhs[32]*(params.A[332])-rhs[33]*(params.A[333])-rhs[34]*(params.A[334])-rhs[35]*(params.A[335])-rhs[36]*(params.A[336])-rhs[37]*(params.A[337])-rhs[38]*(params.A[338])-rhs[39]*(params.A[339])-rhs[40]*(params.A[340])-rhs[41]*(params.A[341])-rhs[42]*(params.A[342])-rhs[43]*(params.A[343])-rhs[44]*(params.A[344])-rhs[45]*(params.A[345])-rhs[46]*(params.A[346])-rhs[47]*(params.A[347])-rhs[48]*(params.A[348])-rhs[49]*(params.A[349])-rhs[50]*(params.A[350])-rhs[51]*(params.A[351])-rhs[52]*(params.A[352])-rhs[53]*(params.A[353])-rhs[54]*(params.A[354])-rhs[55]*(params.A[355])-rhs[56]*(params.A[356])-rhs[57]*(params.A[357])-rhs[58]*(params.A[358])-rhs[59]*(params.A[359]);
  lhs[6] = -rhs[0]*(params.A[360])-rhs[1]*(params.A[361])-rhs[2]*(params.A[362])-rhs[3]*(params.A[363])-rhs[4]*(params.A[364])-rhs[5]*(params.A[365])-rhs[6]*(params.A[366])-rhs[7]*(params.A[367])-rhs[8]*(params.A[368])-rhs[9]*(params.A[369])-rhs[10]*(params.A[370])-rhs[11]*(params.A[371])-rhs[12]*(params.A[372])-rhs[13]*(params.A[373])-rhs[14]*(params.A[374])-rhs[15]*(params.A[375])-rhs[16]*(params.A[376])-rhs[17]*(params.A[377])-rhs[18]*(params.A[378])-rhs[19]*(params.A[379])-rhs[20]*(params.A[380])-rhs[21]*(params.A[381])-rhs[22]*(params.A[382])-rhs[23]*(params.A[383])-rhs[24]*(params.A[384])-rhs[25]*(params.A[385])-rhs[26]*(params.A[386])-rhs[27]*(params.A[387])-rhs[28]*(params.A[388])-rhs[29]*(params.A[389])-rhs[30]*(params.A[390])-rhs[31]*(params.A[391])-rhs[32]*(params.A[392])-rhs[33]*(params.A[393])-rhs[34]*(params.A[394])-rhs[35]*(params.A[395])-rhs[36]*(params.A[396])-rhs[37]*(params.A[397])-rhs[38]*(params.A[398])-rhs[39]*(params.A[399])-rhs[40]*(params.A[400])-rhs[41]*(params.A[401])-rhs[42]*(params.A[402])-rhs[43]*(params.A[403])-rhs[44]*(params.A[404])-rhs[45]*(params.A[405])-rhs[46]*(params.A[406])-rhs[47]*(params.A[407])-rhs[48]*(params.A[408])-rhs[49]*(params.A[409])-rhs[50]*(params.A[410])-rhs[51]*(params.A[411])-rhs[52]*(params.A[412])-rhs[53]*(params.A[413])-rhs[54]*(params.A[414])-rhs[55]*(params.A[415])-rhs[56]*(params.A[416])-rhs[57]*(params.A[417])-rhs[58]*(params.A[418])-rhs[59]*(params.A[419]);
  lhs[7] = -rhs[0]*(params.A[420])-rhs[1]*(params.A[421])-rhs[2]*(params.A[422])-rhs[3]*(params.A[423])-rhs[4]*(params.A[424])-rhs[5]*(params.A[425])-rhs[6]*(params.A[426])-rhs[7]*(params.A[427])-rhs[8]*(params.A[428])-rhs[9]*(params.A[429])-rhs[10]*(params.A[430])-rhs[11]*(params.A[431])-rhs[12]*(params.A[432])-rhs[13]*(params.A[433])-rhs[14]*(params.A[434])-rhs[15]*(params.A[435])-rhs[16]*(params.A[436])-rhs[17]*(params.A[437])-rhs[18]*(params.A[438])-rhs[19]*(params.A[439])-rhs[20]*(params.A[440])-rhs[21]*(params.A[441])-rhs[22]*(params.A[442])-rhs[23]*(params.A[443])-rhs[24]*(params.A[444])-rhs[25]*(params.A[445])-rhs[26]*(params.A[446])-rhs[27]*(params.A[447])-rhs[28]*(params.A[448])-rhs[29]*(params.A[449])-rhs[30]*(params.A[450])-rhs[31]*(params.A[451])-rhs[32]*(params.A[452])-rhs[33]*(params.A[453])-rhs[34]*(params.A[454])-rhs[35]*(params.A[455])-rhs[36]*(params.A[456])-rhs[37]*(params.A[457])-rhs[38]*(params.A[458])-rhs[39]*(params.A[459])-rhs[40]*(params.A[460])-rhs[41]*(params.A[461])-rhs[42]*(params.A[462])-rhs[43]*(params.A[463])-rhs[44]*(params.A[464])-rhs[45]*(params.A[465])-rhs[46]*(params.A[466])-rhs[47]*(params.A[467])-rhs[48]*(params.A[468])-rhs[49]*(params.A[469])-rhs[50]*(params.A[470])-rhs[51]*(params.A[471])-rhs[52]*(params.A[472])-rhs[53]*(params.A[473])-rhs[54]*(params.A[474])-rhs[55]*(params.A[475])-rhs[56]*(params.A[476])-rhs[57]*(params.A[477])-rhs[58]*(params.A[478])-rhs[59]*(params.A[479]);
  lhs[8] = -rhs[0]*(params.A[480])-rhs[1]*(params.A[481])-rhs[2]*(params.A[482])-rhs[3]*(params.A[483])-rhs[4]*(params.A[484])-rhs[5]*(params.A[485])-rhs[6]*(params.A[486])-rhs[7]*(params.A[487])-rhs[8]*(params.A[488])-rhs[9]*(params.A[489])-rhs[10]*(params.A[490])-rhs[11]*(params.A[491])-rhs[12]*(params.A[492])-rhs[13]*(params.A[493])-rhs[14]*(params.A[494])-rhs[15]*(params.A[495])-rhs[16]*(params.A[496])-rhs[17]*(params.A[497])-rhs[18]*(params.A[498])-rhs[19]*(params.A[499])-rhs[20]*(params.A[500])-rhs[21]*(params.A[501])-rhs[22]*(params.A[502])-rhs[23]*(params.A[503])-rhs[24]*(params.A[504])-rhs[25]*(params.A[505])-rhs[26]*(params.A[506])-rhs[27]*(params.A[507])-rhs[28]*(params.A[508])-rhs[29]*(params.A[509])-rhs[30]*(params.A[510])-rhs[31]*(params.A[511])-rhs[32]*(params.A[512])-rhs[33]*(params.A[513])-rhs[34]*(params.A[514])-rhs[35]*(params.A[515])-rhs[36]*(params.A[516])-rhs[37]*(params.A[517])-rhs[38]*(params.A[518])-rhs[39]*(params.A[519])-rhs[40]*(params.A[520])-rhs[41]*(params.A[521])-rhs[42]*(params.A[522])-rhs[43]*(params.A[523])-rhs[44]*(params.A[524])-rhs[45]*(params.A[525])-rhs[46]*(params.A[526])-rhs[47]*(params.A[527])-rhs[48]*(params.A[528])-rhs[49]*(params.A[529])-rhs[50]*(params.A[530])-rhs[51]*(params.A[531])-rhs[52]*(params.A[532])-rhs[53]*(params.A[533])-rhs[54]*(params.A[534])-rhs[55]*(params.A[535])-rhs[56]*(params.A[536])-rhs[57]*(params.A[537])-rhs[58]*(params.A[538])-rhs[59]*(params.A[539]);
  lhs[9] = -rhs[0]*(params.A[540])-rhs[1]*(params.A[541])-rhs[2]*(params.A[542])-rhs[3]*(params.A[543])-rhs[4]*(params.A[544])-rhs[5]*(params.A[545])-rhs[6]*(params.A[546])-rhs[7]*(params.A[547])-rhs[8]*(params.A[548])-rhs[9]*(params.A[549])-rhs[10]*(params.A[550])-rhs[11]*(params.A[551])-rhs[12]*(params.A[552])-rhs[13]*(params.A[553])-rhs[14]*(params.A[554])-rhs[15]*(params.A[555])-rhs[16]*(params.A[556])-rhs[17]*(params.A[557])-rhs[18]*(params.A[558])-rhs[19]*(params.A[559])-rhs[20]*(params.A[560])-rhs[21]*(params.A[561])-rhs[22]*(params.A[562])-rhs[23]*(params.A[563])-rhs[24]*(params.A[564])-rhs[25]*(params.A[565])-rhs[26]*(params.A[566])-rhs[27]*(params.A[567])-rhs[28]*(params.A[568])-rhs[29]*(params.A[569])-rhs[30]*(params.A[570])-rhs[31]*(params.A[571])-rhs[32]*(params.A[572])-rhs[33]*(params.A[573])-rhs[34]*(params.A[574])-rhs[35]*(params.A[575])-rhs[36]*(params.A[576])-rhs[37]*(params.A[577])-rhs[38]*(params.A[578])-rhs[39]*(params.A[579])-rhs[40]*(params.A[580])-rhs[41]*(params.A[581])-rhs[42]*(params.A[582])-rhs[43]*(params.A[583])-rhs[44]*(params.A[584])-rhs[45]*(params.A[585])-rhs[46]*(params.A[586])-rhs[47]*(params.A[587])-rhs[48]*(params.A[588])-rhs[49]*(params.A[589])-rhs[50]*(params.A[590])-rhs[51]*(params.A[591])-rhs[52]*(params.A[592])-rhs[53]*(params.A[593])-rhs[54]*(params.A[594])-rhs[55]*(params.A[595])-rhs[56]*(params.A[596])-rhs[57]*(params.A[597])-rhs[58]*(params.A[598])-rhs[59]*(params.A[599]);
}
void multbyP(double *lhs, double *rhs) {
  /* TODO use the fact that P is symmetric? */
  /* TODO check doubling / half factor etc. */
  lhs[0] = 0;
  lhs[1] = 0;
  lhs[2] = 0;
  lhs[3] = 0;
  lhs[4] = 0;
  lhs[5] = 0;
  lhs[6] = 0;
  lhs[7] = 0;
  lhs[8] = 0;
  lhs[9] = 0;
}
void fillq(void) {
  work.q[0] = 0;
  work.q[1] = 0;
  work.q[2] = 0;
  work.q[3] = 0;
  work.q[4] = 0;
  work.q[5] = 0;
  work.q[6] = 0;
  work.q[7] = 0;
  work.q[8] = 0;
  work.q[9] = 0;
}
void fillh(void) {
  work.h[0] = params.b[0];
  work.h[1] = params.b[1];
  work.h[2] = params.b[2];
  work.h[3] = params.b[3];
  work.h[4] = params.b[4];
  work.h[5] = params.b[5];
  work.h[6] = params.b[6];
  work.h[7] = params.b[7];
  work.h[8] = params.b[8];
  work.h[9] = params.b[9];
  work.h[10] = params.b[10];
  work.h[11] = params.b[11];
  work.h[12] = params.b[12];
  work.h[13] = params.b[13];
  work.h[14] = params.b[14];
  work.h[15] = params.b[15];
  work.h[16] = params.b[16];
  work.h[17] = params.b[17];
  work.h[18] = params.b[18];
  work.h[19] = params.b[19];
  work.h[20] = params.b[20];
  work.h[21] = params.b[21];
  work.h[22] = params.b[22];
  work.h[23] = params.b[23];
  work.h[24] = params.b[24];
  work.h[25] = params.b[25];
  work.h[26] = params.b[26];
  work.h[27] = params.b[27];
  work.h[28] = params.b[28];
  work.h[29] = params.b[29];
  work.h[30] = params.b[30];
  work.h[31] = params.b[31];
  work.h[32] = params.b[32];
  work.h[33] = params.b[33];
  work.h[34] = params.b[34];
  work.h[35] = params.b[35];
  work.h[36] = params.b[36];
  work.h[37] = params.b[37];
  work.h[38] = params.b[38];
  work.h[39] = params.b[39];
  work.h[40] = params.b[40];
  work.h[41] = params.b[41];
  work.h[42] = params.b[42];
  work.h[43] = params.b[43];
  work.h[44] = params.b[44];
  work.h[45] = params.b[45];
  work.h[46] = params.b[46];
  work.h[47] = params.b[47];
  work.h[48] = params.b[48];
  work.h[49] = params.b[49];
  work.h[50] = params.b[50];
  work.h[51] = params.b[51];
  work.h[52] = params.b[52];
  work.h[53] = params.b[53];
  work.h[54] = params.b[54];
  work.h[55] = params.b[55];
  work.h[56] = params.b[56];
  work.h[57] = params.b[57];
  work.h[58] = params.b[58];
  work.h[59] = params.b[59];
}
void fillb(void) {
}
void pre_ops(void) {
}
