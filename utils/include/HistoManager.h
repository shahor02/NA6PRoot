// NA6PCCopyright

// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef _NA6P_HISTOMANAGER_H_
#define _NA6P_HISTOMANAGER_H_

#include <string>
#include <Rtypes.h>
#include "TObjArray.h"

class TH1;
class TH2;
class TH1F;
class TH2F;
class TProfile;
class TGraph;
class TFile;
class TFrame;
class TPad;
class TLegend;
class TVirtualPad;
class TPaveStats;
class TLatex;
class TLegendEntry;

namespace na6p
{

class HistoManager : public TObjArray
{
 public:
  HistoManager(const std::string& dirname = "", const std::string& fname = "histoman.root", bool load = kFALSE, const std::string& prefix = "");
  ~HistoManager() override { Delete(); }

  HistoManager* createClone(const std::string& prefix) const;
  void addPrefix(const std::string& pref = "");

  int getNHistos() const { return mNHistos; }
  TGraph* getGraph(int id) const;
  TH1* getHisto(int id) const;
  TH1* getHisto(const std::string& name) const;
  TH1F* getHisto1F(int id) const;
  TH2F* getHisto2F(int id) const;
  TProfile* getHistoP(int id) const;

  int addHisto(TH1* histo, int at = -1);
  int addGraph(TGraph* gr, int at = -1);
  void delHisto(int at);

  void setFile(TFile* file);
  void setFileName(const std::string& fname);
  const std::string& getFileName() const { return mDefName; }
  void setDirName(const std::string& name) { mDirName = name; }
  const std::string& getDirName() const { return mDirName; }

  void reset();
  void write(TFile* file = nullptr);
  int write(const std::string& flname)
  {
    setFileName(flname);
    write();
    return 0;
  }

  void addHistos(const HistoManager* hm, Double_t c1 = 1.);
  void divideHistos(const HistoManager* hm);
  void multiplyHistos(const HistoManager* hm);
  void scaleHistos(Double_t c1 = 1.);
  void setColor(int tcolor = 1);
  void setMarkerStyle(Style_t mstyle = 1, Size_t msize = 1);
  void setMarkerSize(Size_t msize = 1);
  void sumw2();
  int load(const std::string& fname, const std::string& dirname = "");

  void purify(bool emptyToo = kFALSE);

  void print(Option_t* option = "") const;
  void Clear(Option_t* option = "") override;
  void Delete(Option_t* option = "") override;
  void Compress() override;

  // static methods for histos manipulation
  static TH1* cumulate(TH1* histo, const char* copyName = "_cumul", bool doErr = false);
  static TH1* getBaseHisto(TPad* pad = 0);
  static TH1* getHistosMinMaxRange(TVirtualPad* pad, float& mn, float& mx);
  static void setHistosMinMaxRange(TVirtualPad* pad, float mn, float mx, float marginH = 0.15, float marginL = 0.15);
  static TFrame* getFrame(TPad* pad = 0);
  static TPaveStats* setStatPad(TH1* hst, float x1, float x2, float y1, float y2, Int_t stl = -1, Int_t col = -1);
  static TPaveStats* getStatPad(TH1* hst);
  static void setHStyle(TH1* hst, int col = kRed, int mark = 20, float mrsize = 0.7);
  static void setGStyle(TGraph* hst, int col = kRed, int mark = 20, float mrsize = 0.7);
  static TH1* prof2TH1(TProfile* prof, const char* addName = 0);
  static TH1* invertHisto(TH1* h);
  static TLatex* addLabel(const char* txt, float x = 0.1, float y = 0.9, int color = kBlack, float size = 0.04);
  static TLegendEntry* addLegendEntry(TLegend* leg, const TH1* obj, const char* label = "", Option_t* option = "lpf");
  static TLegendEntry* addLegendEntry(TLegend* leg, const TGraph* obj, const char* label = "", Option_t* option = "lpf");

 private:
  int mNHistos{0};        //! Number of histograms defined
  std::string mDefName{}; //! Default file name
  std::string mDirName{}; //! Directory name in the output file

  ClassDefOverride(HistoManager, 0);
};

} // namespace na6p

#endif // _NA6P_HISTOMANAGER_H_
