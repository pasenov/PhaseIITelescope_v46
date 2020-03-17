//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: HistoManager.cc 72242 2013-07-12 08:44:19Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>
#include <sstream>

#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  :fFileName("phaseiiv1")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);   //enable inactivation of histograms
  
  // Define histograms start values
  //const G4int kNbofStrips = 1016;
  const G4int kMaxHisto = 358;
  const G4int kMaxHisto2 = 37;

  const G4String id[] = { "00", "01", "02", "03" , "04", "05", "06", "07", "08", "09" , "010", "011", "012", "013", "014", "015", "016", "017", "018", "019", "020", "021", "022", "023", "024", "025", "026", "027", "028", "029", "030", "031", "032", "033", "034", "035", "036", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100", "101", "102", "103", "104", "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115", "116", "117", "118", "119", "120", "121", "122", "123", "124", "125", "126", "127", "128", "129", "130", "131", "132", "133", "134", "135", "136", "137", "138", "139", "140", "141", "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172", "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186", "187", "188", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198", "199", "200", "201", "202", "203", "204", "205", "206", "207", "208", "209", "210", "211", "212", "213", "214", "215", "216", "217", "218", "219", "220", "221", "222", "223", "224", "225", "226", "227", "228", "229", "230", "231", "232", "233", "234", "235", "236", "237", "238", "239", "240", "241", "242", "243", "244", "245", "246", "247", "248", "249", "250", "251", "252", "253", "254", "255", "256", "257", "258", "259", "260", "261", "262", "263", "264", "265", "266", "267", "268", "269", "270", "271", "272", "273", "274", "275", "276", "277", "278", "279", "280", "281", "282", "283", "284", "285", "286", "287", "288", "289", "290", "291", "292", "293", "294", "295", "296", "297", "298", "299", "300", "301", "302", "303", "304", "305", "306", "307", "308", "309", "310", "311", "312", "313", "314", "315", "316", "317", "318", "319", "320", "321", "322", "323", "324", "325", "326", "327", "328", "329", "330", "331", "332", "333", "334", "335", "336", "337", "338", "339", "340", "341", "342", "343", "344", "345", "346", "347", "348", "349", "350", "351", "352", "353", "354", "355", "356", "357" };
  const G4String id2[] = { "0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36" };

  const G4String title[] =
                { "dummy",                                                     //0
                  "continuous energy loss along primary track for Sensor 1 (simulation)",   //1
                  "continuous energy loss along primary track for Sensor 2 (simulation)",   //2
                  "energy from secondaries for Sensor 1 (simulation)",                      //3
                  "energy from secondaries for Sensor 2 (simulation)",                      //4
                  "energy from tertiaries for Sensor 1 (simulation)",                       //5
                  "energy from tertiaries for Sensor 2 (simulation)",                       //6
                  "total energy lost by primaries for Sensor 1 (simulation)",           //7
                  "total energy lost by primaries for Sensor 2 (simulation)",           //8
                  "energy spectrum of secondary e-+ for Sensor 1 (simulation)",             //9
                  "energy spectrum of secondary e-+ for Sensor 2 (simulation)",             //10
                  "energy spectrum of secondary gamma for Sensor 1 (simulation)",           //11
                  "energy spectrum of secondary gamma for Sensor 2 (simulation)",           //12
                  "energy spectrum of tertiary e-+ for Sensor 1 (simulation)",              //13
                  "energy spectrum of tertiary e-+ for Sensor 2 (simulation)",              //14
                  "energy spectrum of tertiary gamma for Sensor 1 (simulation)",            //15
                  "energy spectrum of tertiary gamma for Sensor 2 (simulation)",            //16
                  "step size for primary (simulation)",                                       //17
                  "step size for secondaries (simulation)",                                   //18
                  "step size for tertiaries (simulation)",                                    //19
                  "step size for secondaries created in sensors (simulation)",              //20
                  "x-polarization of secondaries (simulation)",                               //21
                  "y-polarization of secondaries (simulation)",                               //22
                  "z-polarization of secondaries (simulation)",                               //23
                  "x-polarization of tertiaries (simulation)",                                //24
                  "y-polarization of tertiaries (simulation)",                                //25
                  "z-polarization of tertiaries (simulation)",                                //26
                  "track length for primary (simulation)",                                    //27
                  "total track length of secondaries (simulation)",                           //28
                  "total track length of tertiaries (simulation)",                            //29
                  "total track length of secondaries created in sensors (simulation)",      //30
                  "track length of secondaries created in sensors (simulation)",            //31
                  "number of positrons created in Sensor 1 (simulation)",       		 //32
                  "number of positrons created in Sensor 2 (simulation)",       		 //33
                  "number of electrons created in Sensor 1 (simulation)",       		 //34
                  "number of electrons created in Sensor 2 (simulation)",       		 //35
                  "number of pi+ created in Sensor 1 (simulation)",             		 //36
                  "number of pi+ created in Sensor 2 (simulation)",             		 //37
                  "number of pi- created in Sensor 1 (simulation)",              		 //38
                  "number of pi- created in Sensor 2 (simulation)",             		 //39
                  "number of mu+ created in Sensor 1 (simulation)",             		 //40
                  "number of mu+ created in Sensor 2 (simulation)",             		 //41
                  "number of mu- created in Sensor 1 (simulation)",             		 //42
                  "number of mu- created in Sensor 2 (simulation)",              		 //43
		  "number of delta e- created in Sensor 1 that have reached Sensor 2 (simulation)", //44
		  "number of delta e- created in Sensor 2 that have reached Sensor 1 (simulation)", //45
                  "number of e+ created before Sensor 1 (simulation)",               	 //46
                  "number of e- created before Sensor 1 (simulation)",               	 //47
		  "number of e+ createdfile:///usr/share/doc/HTML/index.html between Sensor 1 and Sensor 2 (simulation)", 	 //48
		  "number of e- created between Sensor 1 and Sensor 2 (simulation)", 	 //49
                  "number of e+ created after Sensor 2 (simulation)",               	 //50
                  "number of e- created after Sensor 2 (simulation)",                	 //51
                  "continuous energy loss along primary track for strip 505 a Sensor 1 (simulation)",   //52
                  "continuous energy loss along primary track for strip 506 a Sensor 1 (simulation)",   //53
                  "continuous energy loss along primary track for strip 507 a Sensor 1 (simulation)",   //54
                  "continuous energy loss along primary track for strip 508 a Sensor 1 (simulation)",   //55
                  "continuous energy loss along primary track https://root.cern.ch/root/htmldoc/guides/users-guide/Histograms.html#creating-histogramsfor strip 509 a Sensor 1 (simulation)",   //56
                  "continuous energy loss along primary track for strip 510 a Sensor 1 (simulation)",   //57
                  "continuous energy loss along primary track for strip 511 a Sensor 1 (simulation)",   //58
                  "continuous energy loss along primary track for strip 512 a Sensor 1 (simulation)",   //59
                  "continuous energy loss along primary track for strip 505 b Sensor 1 (simulation)",   //60
                  "continuous energy loss along primary track for strip 506 b Sensor 1 (simulation)",   //61
                  "continuous energy loss along primary track for strip 507 b Sensor 1 (simulation)",   //62
                  "continuous energy loss along primary track for strip 508 b Sensor 1 (simulation)",   //63
                  "continuous energy loss along primary track for strip 509 b Sensor 1 (simulation)",   //64
                  "continuous energy loss along primary track for strip 510 b Sensor 1 (simulation)",   //65
                  "continuous energy loss along primary track for strip 511 b Sensor 1 (simulation)",   //66
                  "continuous energy loss along primary track for strip 512 b Sensor 1 (simulation)",   //67
                  "continuous energy loss along primary track for strip 505 a Sensor 2 (simulation)",   //68
                  "continuous energy loss along primary track for strip 506 a Sensor 2 (simulation)",   //69
                  "continuous energy loss along primary track for strip 507 a Sensor 2 (simulation)",   //70
                  "continuous energy loss along primary track for strip 508 a Sensor 2 (simulation)",   //71
                  "continuous energy loss along primary track for strip 509 a Sensor 2 (simulation)",   //72
                  "continuous energy loss along primary track for strip 510 a Sensor 2 (simulation)",   //73
                  "continuous energy loss along primary track for strip 511 a Sensor 2 (simulation)",   //74
                  "continuous energy loss along primary track for strip 512 a Sensor 2 (simulation)",   //75
                  "continuous energy loss along primary track for strip 505 b Sensor 2 (simulation)",   //76
                  "continuous energy loss along primary track for strip 506 b Sensor 2 (simulation)",   //77
                  "continuous energy loss along primary track for strip 507 b Sensor 2 (simulation)",   //78
                  "continuous energy loss along primary track for strip 508 b Sensor 2 (simulation)",   //79
                  "continuous energy loss along primary track for strip 509 b Sensor 2 (simulation)",   //80
                  "continuous energy loss along primary track for strip 510 b Sensor 2 (simulation)",   //81
                  "continuous energy loss along primary track for strip 511 b Sensor 2 (simulation)",   //82
                  "continuous energy loss along primary track for strip 512 b Sensor 2 (simulation)",   //83
                  "number of electrons created in strip 505 a Sensor 1 (simulation)",   //84
                  "number of electrons created in strip 506 a Sensor 1 (simulation)",   //85
                  "number of electrons created in strip 507 a Sensor 1 (simulation)",   //86
                  "number of electrons created in strip 508 a Sensor 1 (simulation)",   //87
                  "number of electrons created in strip 509 a Sensor 1 (simulation)",   //88
                  "number of electrons created in strip 510 a Sensor 1 (simulation)",   //89
                  "number of electrons created in strip 511 a Sensor 1 (simulation)",   //90
                  "number of electrons created in strip 512 a Sensor 1 (simulation)",   //91
                  "number of electrons created in strip 505 b Sensor 1 (simulation)",   //92
                  "number of electrons created in strip 506 b Sensor 1 (simulation)",   //93
                  "number of electrons created in strip 507 b Sensor 1 (simulation)",   //94
                  "number of electrons created in strip 508 b Sensor 1 (simulation)",   //95
                  "number of electrons created in strip 509 b Sensor 1 (simulation)",   //96
                  "number of electrons created in strip 510 b Sensor 1 (simulation)",   //97
                  "number of electrons created in strip 511 b Sensor 1 (simulation)",   //98
                  "number of electrons created in strip 512 b Sensor 1 (simulation)",   //99
                  "number of electrons created in strip 505 a Sensor 2 (simulation)",   //100
                  "number of electrons file:///usr/share/doc/HTML/index.htmlcreated in strip 506 a Sensor 2 (simulation)",   //101
                  "number of electrons created in strip 507 a Sensor 2 (simulation)",   //102
                  "number of electrons created in strip 508 a Sensor 2 (simulation)",   //103
                  "number of electrons created in strip 509 a Sensor 2 (simulation)",   //104
                  "number of electrons created in strip 510 a Sensor 2 (simulation)",   //105
                  "number of electrons created in strip 511 a Sensor 2 (simulation)",   //106
                  "number of electrons created in strip 512 a Sensor 2 (simulation)",   //107
                  "number of electrons created in strip 505 b Sensor 2 (simulation)",   //108
                  "number of electrons created in strip 506 b Sensor 2 (simulation)",   //109
                  "number of electrons created in strip 507 b Sensor 2 (simulation)",   //110
                  "number of electrons created in strip 508 b Sensor 2 (simulation)",   //111
                  "number of electrons created in strip 509 b Sensor 2 (simulation)",   //112
                  "number of electrons created in strip 510 b Sensor 2 (simulation)",   //113
                  "number of electrons created in strip 511 b Sensor 2 (simulation)",   //114
                  "number of electrons created in strip 512 b Sensor 2 (simulation)",   //115
                  "deflection angle of primary particle (rad) (simulation)",     //116
                  "distance between midpoint of primary inside Sensor 1 and AJ (simulation)",   //117
                  "distance between midpoint of primary inside Sensor 2 and AJ (simulation)",   //118
                  "2S first-plane residuals in global x-direction (E'Ex) (simulation)",   //119
                  "2S first-plane residuals in global y-direction (E'Ey) (simulation)",   //120
                  "2S second-plane residuals in global x-direction (F'Fx) (simulation)",   //121
                  "2S second-plane residuals in global y-direction (F'Fy) (simulation)",   //122
                  "number of hits for strip 505 a Sensor 1 (simulation)",   //123
                  "number of hits for strip 506 a Sensor 1 (simulation)",   //124
                  "number of hits for strip 507 a Sensor 1 (simulation)",   //125
                  "number of hits for strip 508 a Sensor 1 (simulation)",   //126
                  "number of hits for strip 509 a Sensor 1 (simulation)",   //127
                  "number of hits for strip 510 a Sensor 1 (simulation)",   //128
                  "number of hits for strip 511 a Sensor 1 (simulation)",   //129
                  "number of hits for strip 512 a Sensor 1 (simulation)",   //130
                  "number of hits for strip 505 b Sensor 1 (simulation)",   //131
                  "number of hits for strip 506 b Sensor 1 (simulation)",   //132
                  "number of hits for strip 507 b Sensor 1 (simulation)",   //133
                  "number of hits for strip 508 b Sensor 1 (simulation)",   //134
                  "number of hits for strip 509 b Sensor 1 (simulation)",   //135
                  "number of hits for strip 510 b Sensor 1 (simulation)",   //136
                  "number of hits for strip 511 b Sensor 1 (simulation)",   //137
                  "number of hits for strip 512 b Sensor 1 (simulation)",   //138
                  "number of hits for strip 505 a Sensor 2 (simulation)",   //139
                  "number of hits for strip 506 a Sensor 2 (simulation)",   //140
                  "number of hits for strip 507 a Sensor 2 (simulation)",   //141
                  "number of hits for strip 508 a Sensor 2 (simulation)",   //142
                  "number of hits for strip 509 a Sensor 2 (simulation)",   //143
                  "number of hits for strip 510 a Sensor 2 (simulation)",   //144
                  "number of hits for strip 511 a Sensor 2 (simulation)",   //145
                  "number of hits for strip 512 a Sensor 2 (simulation)",   //146
                  "number of hits for strip 505 b Sensor 2 (simulation)",   //147
                  "number of hits for strip 506 b Sensor 2 (simulation)",   //148
                  "number of hits for strip 507 b Sensor 2 (simulation)",   //149
                  "number of hits for strip 508 b Sensor 2 (simulation)",   //150
                  "number of hits for strip 509 b Sensor 2 (simulation)",   //151
                  "number of hits for strip 510 b Sensor 2 (simulation)",   //152
                  "number of hits for strip 511 b Sensor 2 (simulation)",   //153
                  "number of hits for strip 512 b Sensor 2 (simulation)",   //154
                  "charge in Sensor 1 corresponding to the deposited energy in Sensor 1 (simulation)",   //155
                  "charge in Sensor 2 corresponding to the deposited energy in Sensor 2 (simulation)",   //156
                  "charge in strip 505 a Sensor 1 (simulation)",   //157
                  "charge in strip 506 a Sensor 1 (simulation)",   //158
                  "charge in strip 507 a Sensor 1 (simulation)",   //159
                  "charge in strip 508 a Sensor 1 (simulation)",   //160
                  "charge in strip 509 a Sensor 1 (simulation)",   //161
                  "charge in strip 510 a Sensor 1 (simulation)",   //162
                  "charge in strip 511 a Sensor 1 (simulation)",   //163
                  "charge in strip 512 a Sensor 1 (simulation)",   //164
                  "charge in strip 505 b Sensor 1 (simulation)",   //165
                  "charge in strip 506 b Sensor 1 (simulation)",   //166
                  "charge in strip 507 b Sensor 1 (simulation)",   //167
                  "charge in strip 508 b Sensor 1 (simulation)",   //168
                  "charge in strip 509 b Sensor 1 (simulation)",   //169
                  "charge in strip 510 b Sensor 1 (simulation)",   //170
                  "charge in strip 511 b Sensor 1 (simulation)",   //171
                  "charge in strip 512 b Sensor 1 (simulation)",   //172
                  "charge in strip 505 a Sensor 2 (simulation)",   //173
                  "charge in strip 506 a Sensor 2 (simulation)",   //174
                  "charge in strip 507 a Sensor 2 (simulation)",   //175
                  "charge in strip 508 a Sensor 2 (simulation)",   //176
                  "charge in strip 509 a Sensor 2 (simulation)",   //177
                  "charge in strip 510 a Sensor 2 (simulation)",   //178
                  "charge in strip 511 a Sensor 2 (simulation)",   //179
                  "charge in strip 512 a Sensor 2 (simulation)",   //180
                  "charge in strip 505 b Sensor 2 (simulation)",   //181
                  "charge in strip 506 b Sensor 2 (simulation)",   //182
                  "charge in strip 507 b Sensor 2 (simulation)",   //183
                  "charge in strip 508 b Sensor 2 (simulation)",   //184
                  "charge in strip 509 b Sensor 2 (simulation)",   //185
                  "charge in strip 510 b Sensor 2 (simulation)",   //186
                  "charge in strip 511 b Sensor 2 (simulation)",   //187
                  "charge in strip 512 b Sensor 2 (simulation)",   //188
                  "cluster size per event: Sensor 1 (simulation)",   //189
                  "cluster size per event: Sensor 2 (simulation)",   //190
                  "reconstructed x-position of primary in sensor 1 (simulation)",    //191
                  "reconstructed y-position of primary in sensor 1 (simulation)",    //192
                  "reconstructed z-position of primary in sensor 1 (simulation)",    //193
                  "reconstructed x-position of primary in sensor 2 (simulation)",    //194
                  "reconstructed y-position of primary in sensor 2 (simulation)",    //195
                  "reconstructed z-position of primary in sensor 2 (simulation)",    //196
                  "cluster size per event: BPIX module 1 (simulation)",    //197
                  "cluster size per event: BPIX module 2 (simulation)",    //198
                  "cluster size per event: BPIX module 3 (simulation)",    //199
                  "cluster size per event: BPIX module 4 (simulation)",    //200
                  "cluster size per event: BPIX module 5 (simulation)",    //201
                  "cluster size per event: BPIX module 6 (simulation)",    //202
                  "cluster size per event: BPIX module 7 (simulation)",    //203
                  "cluster size per event: BPIX module 8 (simulation)",    //204
                  "cluster size per event: BPIX module 9 (simulation)",    //205
                  "cluster size per event: BPIX module 10 (simulation)",    //206
                  "cluster size per event: BPIX module 11 (simulation)",    //207
                  "cluster size per event: BPIX module 12 (simulation)",    //208
                  "cluster size per event: BPIX module 13 (simulation)",    //209
                  "cluster size per event: BPIX module 14 (simulation)",    //210
                  "cluster size per event: BPIX module 15 (simulation)",    //211
                  "cluster size per event: BPIX module 16 (simulation)",    //212
                  "x-value of difference: reconstructed - extrapolated points for the 2S 1st plane (simulation)",     //213
                  "y-value of difference: reconstructed - extrapolated points for the 2S 1st plane (simulation)",     //214
                  "x-value of difference: reconstructed - extrapolated points for the 2S 2nd plane (simulation)",     //215
                  "y-value of difference: reconstructed - extrapolated points for the 2S 2nd plane (simulation)",     //216
                  "reconstructed x-position of primary in pixel module 1 (simulation)",     //217
                  "reconstructed y-position of primary in pixel module 1 (simulation)",     //218
                  "reconstructed z-position of primary in pixel module 1 (simulation)",     //219
                  "reconstructed x-position of primary in pixel module 2 (simulation)",     //220
                  "reconstructed y-position of primary in pixel module 2 (simulation)",     //221
                  "reconstructed z-position of primary in pixel module 2 (simulation)",     //222
                  "reconstructed x-position of primary in pixel module 3 (simulation)",     //223
                  "reconstructed y-position of primary in pixel module 3 (simulation)",     //224
                  "reconstructed z-position of primary in pixel module 3 (simulation)",     //225
                  "reconstructed x-position of primary in pixel module 4 (simulation)",     //226
                  "reconstructed y-position of primary in pixel module 4 (simulation)",     //227
                  "reconstructed z-position of primary in pixel module 4 (simulation)",     //228
                  "reconstructed x-position of primary in pixel module 5 (simulation)",     //229
                  "reconstructed y-position of primary in pixel module 5 (simulation)",     //230
                  "reconstructed z-position of primary in pixel module 5 (simulation)",     //231
                  "reconstructed x-position of primary in pixel module 6 (simulation)",     //232
                  "reconstructed y-position of primary in pixel module 6 (simulation)",     //233
                  "reconstructed z-position of primary in pixel module 6 (simulation)",     //234
                  "reconstructed x-position of primary in pixel module 7 (simulation)",     //235
                  "reconstructed y-position of primary in pixel module 7 (simulation)",     //236
                  "reconstructed z-position of primary in pixel module 7 (simulation)",     //237
                  "reconstructed x-position of primary in pixel module 8 (simulation)",     //238
                  "reconstructed y-position of primary in pixel module 8 (simulation)",     //239
                  "reconstructed z-position of primary in pixel module 8 (simulation)",     //240
                  "reconstructed x-position of primary in pixel module 9 (simulation)",     //241
                  "reconstructed y-position of primary in pixel module 9 (simulation)",     //242
                  "reconstructed z-position of primary in pixel module 9 (simulation)",     //243
                  "reconstructed x-position of primary in pixel module 10 (simulation)",    //244
                  "reconstructed y-position of primary in pixel module 10 (simulation)",    //245
                  "reconstructed z-position of primary in pixel module 10 (simulation)",    //246
                  "reconstructed x-position of primary in pixel module 11 (simulation)",    //247
                  "reconstructed y-position of primary in pixel module 11 (simulation)",    //248
                  "reconstructed z-position of primary in pixel module 11 (simulation)",    //249
                  "reconstructed x-position of primary in pixel module 12 (simulation)",    //250
                  "reconstructed y-position of primary in pixel module 12 (simulation)",    //251
                  "reconstructed z-position of primary in pixel module 12 (simulation)",    //252
                  "reconstructed x-position of primary in pixel module 13 (simulation)",    //253
                  "reconstructed y-position of primary in pixel module 13 (simulation)",    //254
                  "reconstructed z-position of primary in pixel module 13 (simulation)",    //255
                  "reconstructed x-position of primary in pixel module 14 (simulation)",    //256
                  "reconstructed y-position of primary in pixel module 14 (simulation)",    //257
                  "reconstructed z-position of primary in pixel module 14 (simulation)",    //258
                  "reconstructed x-position of primary in pixel module 15 (simulation)",    //259
                  "reconstructed y-position of primary in pixel module 15 (simulation)",    //260
                  "reconstructed z-position of primary in pixel module 15 (simulation)",    //261
                  "reconstructed x-position of primary in pixel module 16 (simulation)",    //262
                  "reconstructed y-position of primary in pixel module 16 (simulation)",    //263
                  "reconstructed z-position of primary in pixel module 16 (simulation)",    //264
                  "Telescope first-layer residuals in global x-direction (A'Ax) (simulation)",   //265
                  "Telescope first-layer residuals in global y-direction (A'Ay) (simulation)",   //266
                  "Telescope second-layer residuals in global x-direction (B'Bx) (simulation)",   //267
                  "Telescope second-layer residuals in global y-direction (B'By) (simulation)",   //268
                  "Telescope third-layer residuals in global x-direction (C'Cx) (simulation)",   //269
                  "Telescope third-layer residuals in global y-direction (C'Cy) (simulation)",   //270
                  "Telescope fourth-layer residuals in global x-direction (D'Dx) (simulation)",   //271
                  "Telescope fourth-layer residuals in global y-direction (D'Dy) (simulation)",   //272
                  "Telescope fifth-layer residuals in global x-direction (G'Gx) (simulation)",   //273
                  "Telescope fifth-layer residuals in global y-direction (G'Gy) (simulation)",   //274
                  "Telescope sixth-layer residuals in global x-direction (H'Hx) (simulation)",   //275
                  "Telescope sixth-layer residuals in global y-direction (H'Hy) (simulation)",   //276
                  "Telescope seventh-layer residuals in global x-direction (I'Ix) (simulation)",   //277
                  "Telescope seventh-layer residuals in global y-direction (I'Iy) (simulation)",   //278
                  "Telescope eighth-layer residuals in global x-direction (J'Jx) (simulation)",   //279
                  "Telescope eighth-layer residuals in global y-direction (J'Jy) (simulation)",   //280
                  "2S first-plane residuals in local x-direction (E'Ex) (simulation)",   //281
                  "2S first-plane residuals in local y-direction (E'Ey) (simulation)",   //282
                  "2S second-plane residuals in local x-direction (F'Fx) (simulation)",   //283
                  "2S second-plane residuals in local y-direction (F'Fy) (simulation)",   //284
                  "Telescope first-layer residuals in local x-direction (A'Ax) (simulation)",   //285
                  "Telescope first-layer residuals in local y-direction (A'Ay) (simulation)",   //286
                  "Telescope second-layer residuals in local x-direction (B'Bx) (simulation)",   //287
                  "Telescope second-layer residuals in local y-direction (B'By) (simulation)",   //288
                  "Telescope third-layer residuals in local x-direction (C'Cx) (simulation)",   //289
                  "Telescope third-layer residuals in local y-direction (C'Cy) (simulation)",   //290
                  "Telescope fourth-layer residuals in local x-direction (D'Dx) (simulation)",   //291
                  "Telescope fourth-layer residuals in local y-direction (D'Dy) (simulation)",   //292
                  "Telescope fifth-layer residuals in local x-direction (G'Gx) (simulation)",   //293
                  "Telescope fifth-layer residuals in local y-direction (G'Gy) (simulation)",   //294
                  "Telescope sixth-layer residuals in local x-direction (H'Hx) (simulation)",   //295
                  "Telescope sixth-layer residuals in local y-direction (H'Hy) (simulation)",   //296
                  "Telescope seventh-layer residuals in local x-direction (I'Ix) (simulation)",   //297
                  "Telescope seventh-layer residuals in local y-direction (I'Iy) (simulation)",   //298
                  "Telescope eighth-layer residuals in local x-direction (J'Jx) (simulation)",   //299
                  "Telescope eighth-layer residuals in local y-direction (J'Jy) (simulation)",   //300
                  "continuous energy loss along primary track for telescope layer 1 (simulation)",   //301
                  "energy from secondaries for telescope layer 1 (simulation)",                      //302
                  "total energy lost by primary beam for telescope layer 1 (simulation)",           //303
                  "continuous energy loss along primary track for telescope layer 2 (simulation)",   //304
                  "energy from secondaries for telescope layer 2 (simulation)",                      //305
                  "total energy lost by primary beam for telescope layer 2 (simulation)",           //306
                  "continuous energy loss along primary track for telescope layer 3 (simulation)",   //307
                  "energy from secondaries for telescope layer 3 (simulation)",                      //308
                  "total energy lost by primary beam for telescope layer 3 (simulation)",           //309
                  "continuous energy loss along primary track for telescope layer 4 (simulation)",   //310
                  "energy from secondaries for telescope layer 4 (simulation)",                      //311
                  "total energy lost by primary beam for telescope layer 4 (simulation)",           //312
                  "continuous energy loss along primary track for telescope layer 5 (simulation)",   //313
                  "energy from secondaries for telescope layer 5 (simulation)",                      //314
                  "total energy lost by primary beam for telescope layer 5 (simulation)",           //315
                  "continuous energy loss along primary track for telescope layer 6 (simulation)",   //316
                  "energy from secondaries for telescope layer 6 (simulation)",                      //317
                  "total energy lost by primary beam for telescope layer 6 (simulation)",           //318
                  "continuous energy loss along primary track for telescope layer 7 (simulation)",   //319
                  "energy from secondaries for telescope layer 7 (simulation)",                      //320
                  "total energy lost by primary beam for telescope layer 7 (simulation)",           //321
                  "continuous energy loss along primary track for telescope layer 8 (simulation)",   //322
                  "energy from secondaries for telescope layer 8 (simulation)",                      //323
                  "total energy lost by primary beam for telescope layer 8 (simulation)",           //324
                  "cluster size-X per event: BPIX module 1 (simulation)",    //325
                  "cluster size-X per event: BPIX module 2 (simulation)",    //326
                  "cluster size-X per event: BPIX module 3 (simulation)",    //327
                  "cluster size-X per event: BPIX module 4 (simulation)",    //328
                  "cluster size-X per event: BPIX module 5 (simulation)",    //329
                  "cluster size-X per event: BPIX module 6 (simulation)",    //330
                  "cluster size-X per event: BPIX module 7 (simulation)",    //331
                  "cluster size-X per event: BPIX module 8 (simulation)",    //332
                  "cluster size-X per event: BPIX module 9 (simulation)",    //333
                  "cluster size-X per event: BPIX module 10 (simulation)",    //334
                  "cluster size-X per event: BPIX module 11 (simulation)",    //335
                  "cluster size-X per event: BPIX module 12 (simulation)",    //336
                  "cluster size-X per event: BPIX module 13 (simulation)",    //337
                  "cluster size-X per event: BPIX module 14 (simulation)",    //338
                  "cluster size-X per event: BPIX module 15 (simulation)",    //339
                  "cluster size-X per event: BPIX module 16 (simulation)",    //340
                  "cluster size-Y per event: BPIX module 1 (simulation)",    //341
                  "cluster size-Y per event: BPIX module 2 (simulation)",    //342
                  "cluster size-Y per event: BPIX module 3 (simulation)",    //343
                  "cluster size-Y per event: BPIX module 4 (simulation)",    //344
                  "cluster size-Y per event: BPIX module 5 (simulation)",    //345
                  "cluster size-Y per event: BPIX module 6 (simulation)",    //346
                  "cluster size-Y per event: BPIX module 7 (simulation)",    //347
                  "cluster size-Y per event: BPIX module 8 (simulation)",    //348
                  "cluster size-Y per event: BPIX module 9 (simulation)",    //349
                  "cluster size-Y per event: BPIX module 10 (simulation)",    //350
                  "cluster size-Y per event: BPIX module 11 (simulation)",    //351
                  "cluster size-Y per event: BPIX module 12 (simulation)",    //352
                  "cluster size-Y per event: BPIX module 13 (simulation)",    //353
                  "cluster size-Y per event: BPIX module 14 (simulation)",    //354
                  "cluster size-Y per event: BPIX module 15 (simulation)",    //355
                  "cluster size-Y per event: BPIX module 16 (simulation)",    //356
                  "Cluster charge for M3173 (left module of Layer 1)"     //357
                 };

  const G4String title2[] =
                { "dummy",                                                     //0
                  "Y vs X in global frame for telescope Layer 1 (simulation)",   //1
                  "Y vs X in global frame for telescope Layer 2 (simulation)",   //2
                  "Y vs X in global frame for telescope Layer 3 (simulation)",   //3
                  "Y vs X in global frame for telescope Layer 4 (simulation)",   //4
                  "Y vs X in global frame for telescope Layer 5 (simulation)",   //5
                  "Y vs X in global frame for telescope Layer 6 (simulation)",   //6
                  "Y vs X in global frame for telescope Layer 7 (simulation)",   //7
                  "Y vs X in global frame for telescope Layer 8 (simulation)",   //8
                  "Y vs X in global frame for DUT Sensor 1 (simulation)",   //9
                  "Y vs X in global frame for DUT Sensor 2 (simulation)",   //10
                  "Y_{local} vs X_{local} in local frame for telescope Layer 1 (simulation)",   //11
                  "Y_{local} vs X_{local} in local frame for telescope Layer 2 (simulation)",   //12
                  "Y_{local} vs X_{local} in local frame for telescope Layer 3 (simulation)",   //13
                  "Y_{local} vs X_{local} in local frame for telescope Layer 4 (simulation)",   //14
                  "Y_{local} vs X_{local} in local frame for telescope Layer 5 (simulation)",   //15
                  "Y_{local} vs X_{local} in local frame for telescope Layer 6 (simulation)",   //16
                  "Y_{local} vs X_{local} in local frame for telescope Layer 7 (simulation)",   //17
                  "Y_{local} vs X_{local} in local frame for telescope Layer 8 (simulation)",   //18
                  "Y_{local} vs X_{local} in local frame for DUT Sensor 1 (simulation)",   //19
                  "Y_{local} vs X_{local} in local frame for DUT Sensor 2 (simulation)",   //20
		  "Cluster Occupancy for Module 1 (right module of Layer 1, simulation)",  //21
		  "Cluster Occupancy for Module 2 (left module of Layer 1, simulation)",   //22
		  "Cluster Occupancy for Module 3 (right module of Layer 2, simulation)",  //23
		  "Cluster Occupancy for Module 4 (left module of Layer 2, simulation)",   //24
		  "Cluster Occupancy for Module 5 (right module of Layer 3, simulation)",  //25
		  "Cluster Occupancy for Module 6 (left module of Layer 3, simulation)",   //26
		  "Cluster Occupancy for Module 7 (right module of Layer 4, simulation)",  //27
		  "Cluster Occupancy for Module 8 (left module of Layer 4, simulation)",   //28
		  "Cluster Occupancy for Module 9 (right module of Layer 5, simulation)",  //29
		  "Cluster Occupancy for Module 10 (left module of Layer 5, simulation)",   //30
		  "Cluster Occupancy for Module 11 (right module of Layer 6, simulation)",  //31
		  "Cluster Occupancy for Module 12 (left module of Layer 6, simulation)",   //31
		  "Cluster Occupancy for Module 13 (right module of Layer 7, simulation)",  //33
		  "Cluster Occupancy for Module 14 (left module of Layer 7, simulation)",   //34
		  "Cluster Occupancy for Module 15 (right module of Layer 8, simulation)",  //35
		  "Cluster Occupancy for Module 16 (left module of Layer 8, simulation)"    //36
                 };
            
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 101;
  G4double vmin = -0.5;
  G4double vmax = 100.5;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
  for (G4int k2=0; k2<kMaxHisto2; k2++) {
    G4int ih2 = analysisManager->CreateH2(id2[k2], title2[k2], nbins, vmin, vmax, nbins, vmin, vmax);
    analysisManager->SetH2Activation(ih2, false);
  }

  //G4int ih0 = analysisManager->CreateH2("0", "dummy", nbins, vmin, vmax,  //real point
						      //nbins, vmin, vmax);   //reconstructed point
  //analysisManager->SetH2Activation(ih0, false);

  // nTuples, column l: 1 to 1016 correspond to Sensor 1, 2 to 2032 correspond to Sensor 2
  //
  /*analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetFirstNtupleId(1);       
  analysisManager->CreateNtuple("1", "Primary Particle Tuple 1");
  for (G4int l=0; l<1017; l++) {
     std::ostringstream os1;
     os1 <<"EdepStripSensor1a_" << l;
     analysisManager->CreateNtupleDColumn(os1.str());          //
  }
  for (G4int m1=0; m1<1017; m1++) {
     std::ostringstream os2;
     os2 <<"EdepStripSensor1b_" << m1;
     analysisManager->CreateNtupleDColumn(os2.str());          //
  }
  for (G4int ma=0; ma<1017; ma++) {
     std::ostringstream os3;
     os3 <<"EdepStripSensor2a_" << ma;
     analysisManager->CreateNtupleDColumn(os3.str());          //
  }
  for (G4int mb=0; mb<1017; mb++) {
     std::ostringstream os4;
     os4 <<"EdepStripSensor2b_" << mb;
     analysisManager->CreateNtupleDColumn(os4.str());          //
  }

  for (G4int l1=0; l1<1017; l1++) {
     std::ostringstream os5;
     os5 <<"NbHitsSensor1a_" << l1;
     analysisManager->CreateNtupleDColumn(os5.str());          //
  }
  for (G4int l2=0; l2<1017; l2++) {
     std::ostringstream os6;
     os6 <<"NbHitsSensor1b_" << l2;
     analysisManager->CreateNtupleDColumn(os6.str());          //
  }
  for (G4int l3=0; l3<1017; l3++) {
     std::ostringstream os7;
     os7 <<"NbHitsSensor2a_" << l3;
     analysisManager->CreateNtupleDColumn(os7.str());          //
  }
  for (G4int l4=0; l4<1017; l4++) {
     std::ostringstream os8;
     os8 <<"NbHitsSensor2b_" << l4;
     analysisManager->CreateNtupleDColumn(os8.str());          //
  }

  G4int NbBPIX = 16;
  G4int NbROC = 16;
  G4int NbRows = 52;
  G4int NbCols = 80;

  for (G4int ne1=0; ne1<=NbBPIX; ne1++) {
     for (G4int ne2=0; ne2<=NbROC; ne2++) {
        for (G4int ne3=0; ne3<=NbRows; ne3++) {
           for (G4int ne4=0; ne4<=NbCols; ne4++) {
              std::ostringstream os9;
     	      os9 <<"EdepPixel_" << ne1 <<"_" << ne2 <<"_" << ne3 << "_" << ne4;
     	      //analysisManager->CreateNtupleDColumn(os9.str());  //Name convention: EdepPixel_NbBPIX_NbROC_NbRow_NbCol, zeros are dummy
           }
        }
     }
  }

  for (G4int n1=0; n1<=NbBPIX; n1++) {
     for (G4int n2=0; n2<=NbROC; n2++) {
        for (G4int n3=0; n3<=NbRows; n3++) {
           for (G4int n4=0; n4<=NbCols; n4++) {
              std::ostringstream os10;
     	      os10 <<"NbHitsPixel_" << n1 <<"_" << n2 <<"_" << n3 << "_" << n4;
     	      //analysisManager->CreateNtupleDColumn(os10.str());  //Name convention: NbHitsPixel_NbBPIX_NbROC_NbRow_NbCol, zeros are dummy
           }
        }
     }
  }

  //for (G4int n=0; n<1017; n++) {
     //std::ostringstream os3;
     //os3 <<"NbElectronsStripSensor1_" << n;
     //analysisManager->CreateNtupleDColumn(os3.str());          //
  //}
  //for (G4int o=0; o<1017; o++) {
     //std::ostringstream os4;
     //os4 <<"NbElectronsStripSensor2_" << o;
     //analysisManager->CreateNtupleDColumn(os4.str());          //
  //}
  analysisManager->FinishNtuple();
  
  //analysisManager->SetNtupleActivation(false);  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
