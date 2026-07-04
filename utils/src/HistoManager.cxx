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

#include <TDirectory.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TNamed.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <TPad.h>
#include <TFrame.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TVirtualPad.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "fairlogger/Logger.h"
#include "HistoManager.h"

namespace na6p
{

HistoManager::HistoManager(const std::string& dirname, const std::string& fname, bool load, const std::string& prefix) : mDirName(dirname)
{
  setFileName(fname);
  if (load && !mDefName.empty()) {
    int nh = this->load(fname, dirname);
    LOGP(info, "HistoManager::load was requested: got {} histos from {}/{}", nh, fname, dirname);
    if (!prefix.empty()) {
      addPrefix(prefix);
    }
  }
}

HistoManager* HistoManager::createClone(const std::string& prefix) const
{
  auto* hm = static_cast<HistoManager*>(Clone());
  hm->addPrefix(prefix);
  for (int i = 0; i < GetLast() + 1; ++i) {
    TObject* obj = hm->UncheckedAt(i);
    if (!obj) {
      continue;
    }
    if (auto* histo = dynamic_cast<TH1*>(obj)) {
      histo->SetDirectory(nullptr);
    }
  }
  hm->mNHistos = mNHistos;
  hm->setFileName(mDefName);
  hm->setDirName(mDirName);
  return hm;
}

int HistoManager::addHisto(TH1* histo, int at)
{
  if (!histo) {
    return mNHistos;
  }
  if (at < 0) {
    at = mNHistos;
  }
  AddAtAndExpand(histo, at);
  histo->SetDirectory(nullptr);
  histo->SetUniqueID(at + 1);
  return mNHistos++;
}

TGraph* HistoManager::getGraph(int id) const
{
  return id <= GetLast() ? dynamic_cast<TGraph*>(UncheckedAt(id)) : nullptr;
}

TH1* HistoManager::getHisto(int id) const
{
  return id <= GetLast() ? dynamic_cast<TH1*>(UncheckedAt(id)) : nullptr;
}

TH1* HistoManager::getHisto(const std::string& name) const
{
  return dynamic_cast<TH1*>(FindObject(name.c_str()));
}

TH1F* HistoManager::getHisto1F(int id) const
{
  return dynamic_cast<TH1F*>(UncheckedAt(id));
}

TH2F* HistoManager::getHisto2F(int id) const
{
  return dynamic_cast<TH2F*>(UncheckedAt(id));
}

TProfile* HistoManager::getHistoP(int id) const
{
  return dynamic_cast<TProfile*>(UncheckedAt(id));
}

int HistoManager::addGraph(TGraph* gr, int at)
{
  if (!gr) {
    return mNHistos;
  }
  if (at < 0) {
    at = mNHistos;
  }
  AddAtAndExpand(gr, at);
  gr->SetUniqueID(at + 1);
  return mNHistos++;
}

void HistoManager::Compress()
{
  TObjArray::Compress();
  for (int i = 0; i < GetLast() + 1; ++i) {
    if (TObject* histo = At(i)) {
      histo->SetUniqueID(i + 1);
    }
  }
}

void HistoManager::write(TFile* file)
{
  if (!mNHistos) {
    return;
  }

  bool localFile = kFALSE;
  TFile* lfile = nullptr;
  const char* dirName = nullptr;
  if (file) {
    lfile = file;
  } else {
    auto* tmpF = static_cast<TFile*>(gROOT->GetListOfFiles()->FindObject(mDefName.c_str()));
    if (tmpF && tmpF->IsOpen()) {
      TString opt = tmpF->GetOption();
      opt.ToLower();
      if (!opt.Contains("read")) {
        lfile = tmpF;
        tmpF->cd();
      }
    }
  }

  TString pwd = gDirectory->GetPath();
  if (!lfile) {
    std::string originalName = mDefName;
    if (mDefName.empty() || mDefName[0] == ' ') {
      mDefName = "histoman";
    }
    TString rootName = mDefName.c_str();
    if (!rootName.Contains(".root")) {
      mDefName += ".root";
    }
    lfile = TFile::Open(mDefName.c_str(), "UPDATE");
    mDefName = originalName;
    localFile = kTRUE;
  }

  lfile->cd();
  dirName = mDirName.c_str();
  if (dirName && dirName[0] && dirName[0] != ' ') {
    if (!lfile->Get(dirName)) {
      lfile->mkdir(dirName);
    }
    lfile->cd(dirName);
  }
  LOGP(info, "Writing histograms to: {}/{}", lfile->GetPath(), dirName);

  for (int i = 0; i < GetLast() + 1; ++i) {
    TObject* obj = UncheckedAt(i);
    if (!obj) {
      continue;
    }
    auto* histo = dynamic_cast<TH1*>(obj);
    TDirectory* dr = nullptr;
    if (histo) {
      dr = histo->GetDirectory();
      histo->SetDirectory(nullptr);
    }
    obj->Write(nullptr, TObject::kOverwrite);
    if (dr && histo) {
      histo->SetDirectory(dr);
    }
  }

  if (localFile) {
    lfile->Close();
    delete lfile;
  }
  auto* oldDir = static_cast<TDirectory*>(gROOT->GetListOfFiles()->FindObject(pwd.Data()));
  if (oldDir) {
    oldDir->cd();
  }
}

void HistoManager::Clear(Option_t*)
{
  int nent = GetLast() + 1;
  for (int i = 0; i < nent; ++i) {
    TObject* hh = UncheckedAt(i);
    if (!hh) {
      continue;
    }
    RemoveAt(i);
    --mNHistos;
  }
}

//_______________________________________________________________
void HistoManager::Delete(Option_t*)
{
  int nent = GetLast() + 1;
  for (int i = 0; i < nent; ++i) {
    TObject* hh = UncheckedAt(i);
    if (!hh) {
      continue;
    }
    RemoveAt(i);
    delete hh;
  }
  mNHistos = 0;
}

void HistoManager::print(Option_t* option) const
{
  int nent = GetLast() + 1;
  for (int i = 0; i < nent; ++i) {
    TObject* hh = UncheckedAt(i);
    if (!hh) {
      continue;
    }
    LOGP(info, "At position #{}", i);
    hh->Print(option);
  }
  LOGP(info, "Total number of defined histograms: %d", mNHistos);
  LOGP(info, "Current output path: {}/{}", mDefName, mDirName);
}

void HistoManager::addPrefix(const std::string& pref)
{
  if (pref.empty()) {
    return;
  }
  int nent = GetLast() + 1;
  for (int i = 0; i < nent; ++i) {
    TObject* hh = UncheckedAt(i);
    if (!hh) {
      continue;
    }
    if (hh->InheritsFrom("TNamed")) {
      auto name = pref + hh->GetName();
      static_cast<TNamed*>(hh)->SetName(name.c_str());
    }
  }
}

void HistoManager::addHistos(const HistoManager* hm, Double_t c1)
{
  if (!hm) {
    return;
  }
  int nent = GetLast() + 1;
  int nent1 = hm->GetLast() + 1;
  if (nent != nent1) {
    Error("addHistos", "HistoManagers have different content: %d vs %d", nent, nent1);
    return;
  }
  for (int i = 0; i < nent; ++i) {
    TH1* hh1 = getHisto(i);
    TH1* hh2 = hm->getHisto(i);
    if (!hh1 || !hh2) {
      continue;
    }
    hh1->Add(hh2, c1);
  }
}

void HistoManager::divideHistos(const HistoManager* hm)
{
  if (!hm) {
    return;
  }
  int nent = GetLast() + 1;
  int nent1 = hm->GetLast() + 1;
  if (nent != nent1) {
    Error("divideHistos", "HistoManagers have different content: %d vs %d", nent, nent1);
    return;
  }
  for (int i = 0; i < nent; ++i) {
    TH1* hh1 = getHisto(i);
    TH1* hh2 = hm->getHisto(i);
    if (!hh1 || !hh2) {
      continue;
    }
    hh1->Divide(hh2);
  }
}

//_______________________________________________________________
void HistoManager::multiplyHistos(const HistoManager* hm)
{
  if (!hm) {
    return;
  }
  int nent = GetLast() + 1;
  int nent1 = hm->GetLast() + 1;
  if (nent != nent1) {
    Error("multiplyHistos", "HistoManagers have different content: %d vs %d", nent, nent1);
    return;
  }
  for (int i = 0; i < nent; ++i) {
    TH1* hh1 = getHisto(i);
    TH1* hh2 = hm->getHisto(i);
    if (!hh1 || !hh2) {
      continue;
    }
    hh1->Multiply(hh2);
  }
}

void HistoManager::scaleHistos(Double_t c1)
{
  int nent = GetLast() + 1;
  for (int i = 0; i < nent; ++i) {
    TH1* hh1 = getHisto(i);
    if (hh1) {
      hh1->Scale(c1);
    }
  }
}

void HistoManager::sumw2()
{
  int nent = GetLast() + 1;
  for (int i = 0; i < nent; ++i) {
    auto* hh1 = dynamic_cast<TH1*>(UncheckedAt(i));
    if (hh1) {
      hh1->Sumw2();
    }
  }
}

void HistoManager::setFile(TFile* file)
{
  if (file) {
    mDefName = file->GetName();
  }
}

void HistoManager::delHisto(int at)
{
  TH1* hist = getHisto(at);
  if (hist) {
    RemoveAt(at);
    delete hist;
    --mNHistos;
  }
}

void HistoManager::purify(bool emptyToo)
{
  int last = GetLast() + 1;
  if (emptyToo) {
    for (int i = 0; i < last; ++i) {
      TH1* hist = getHisto(i);
      if (!hist) {
        continue;
      }
      if (hist->GetEntries() < 1) {
        delHisto(i);
      }
    }
  }
  Compress();
}

void HistoManager::setFileName(const std::string& name)
{
  TString sName = name;
  gSystem->ExpandPathName(sName);
  mDefName = sName.Data();
}

void HistoManager::reset()
{
  int last = GetLast() + 1;
  for (int i = 0; i < last; ++i) {
    TH1* hist = getHisto(i);
    if (!hist) {
      continue;
    }
    hist->Reset();
  }
}

int HistoManager::load(const std::string& fname, const std::string& dirname)
{
  TString sName = fname;
  if (gSystem->ExpandPathName(sName)) {
    LOGP(error, "Cannot expand file name {}", fname);
    return 0;
  }
  TFile* file = TFile::Open(sName);
  if (!file) {
    LOGP(error, "No file {}", fname);
    return 0;
  }
  if (!dirname.empty() && dirname[0] != ' ') {
    if (!file->Get(dirname.c_str())) {
      LOGP(error, "No {} directory in file {}", dirname, fname);
      file->Close();
      delete file;
      return 0;
    }
    file->cd(dirname.c_str());
  }
  int count = 0;
  TList* lst = gDirectory->GetListOfKeys();
  TIter nextKey(lst);
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(nextKey()))) {
    if (FindObject(key->GetName())) {
      continue;
    }
    TString clName = key->GetClassName();
    if (!(clName.BeginsWith("TH") || clName.BeginsWith("TProfile") || clName.BeginsWith("TGraph"))) {
      printf("Object %s of type %s is not processed\n", key->GetName(), clName.Data());
      continue;
    }
    TObject* hst = key->ReadObj();
    int id = hst->GetUniqueID();
    if (auto* h = dynamic_cast<TH1*>(hst)) {
      addHisto(h, id - 1);
      ++count;
      continue;
    }
    if (auto* gr = dynamic_cast<TGraph*>(hst)) {
      addGraph(gr, id - 1);
      ++count;
    }
  }
  file->Close();
  delete file;
  auto nm = fname;
  if (!dirname.empty()) {
    nm += fmt::format("/{}", dirname);
  }
  SetName(nm.c_str());
  return count;
}

void HistoManager::setColor(int tcolor)
{
  int last = GetLast() + 1;
  for (int i = 0; i < last; ++i) {
    TH1* hist = getHisto(i);
    if (!hist) {
      continue;
    }
    hist->SetLineColor(tcolor);
    hist->SetMarkerColor(tcolor);
    TList* lst = hist->GetListOfFunctions();
    if (lst) {
      int nf = lst->GetSize();
      for (int j = 0; j < nf; ++j) {
        TObject* fnc = lst->At(j);
        if (fnc->InheritsFrom("TF1")) {
          static_cast<TF1*>(fnc)->SetLineColor(tcolor);
          static_cast<TF1*>(fnc)->SetLineWidth(1);
          static_cast<TF1*>(fnc)->ResetBit(TF1::kNotDraw);
        } else if (fnc->InheritsFrom("TPaveStats")) {
          static_cast<TPaveStats*>(fnc)->SetTextColor(tcolor);
        }
      }
    }
  }
}

void HistoManager::setMarkerStyle(Style_t mstyle, Size_t msize)
{
  int last = GetLast() + 1;
  for (int i = 0; i < last; ++i) {
    TH1* hist = getHisto(i);
    if (!hist) {
      continue;
    }
    hist->SetMarkerStyle(mstyle);
    hist->SetMarkerSize(msize);
  }
}

void HistoManager::setMarkerSize(Size_t msize)
{
  int last = GetLast() + 1;
  for (int i = 0; i < last; ++i) {
    TH1* hist = getHisto(i);
    if (!hist) {
      continue;
    }
    hist->SetMarkerSize(msize);
  }
}

TH1* HistoManager::cumulate(TH1* histo, const char* copyName, bool doErr)
{
  // create cumulative histo
  TString nname = copyName;
  if (nname.IsNull()) {
    nname = histo->GetName();
    nname += "_cml";
  }
  TH1* cml = (TH1*)histo->Clone(nname.Data());
  int nb = histo->GetNbinsX();
  double sm = 0;
  double sme = 0;
  //
  for (int i = 1; i <= nb; i++) {
    sm += histo->GetBinContent(i);
    cml->SetBinContent(i, sm);
    if (!doErr)
      continue;
    double ee = histo->GetBinError(i);
    sme += ee * ee;
    cml->SetBinError(i, sme > 0 ? TMath::Sqrt(sme) : 0.);
  }
  return cml;
}

TH1* HistoManager::getBaseHisto(TPad* pad)
{
  if (!pad)
    pad = (TPad*)gPad;
  if (!pad)
    return 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TH1* hst = 0;
  for (int i = 0; i < size; i++) {
    TObject* obj = lst->At(i);
    if (!obj)
      continue;
    if (obj->InheritsFrom("TH1")) {
      hst = (TH1*)obj;
      break;
    }
  }
  return hst;
}

TH1* HistoManager::getHistosMinMaxRange(TVirtualPad* pad, float& mn, float& mx)
{
  mn = 1e9;
  mx = -1e9;
  if (!pad)
    pad = (TPad*)gPad;
  if (!pad)
    return 0;
  TH1* hbase = 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TH1* hst = 0;
  for (int i = 0; i < size; i++) {
    TObject* obj = lst->At(i);
    if (!obj)
      continue;
    if (obj->InheritsFrom("TH1")) {
      hst = (TH1*)obj;
      float mnh = hst->GetBinContent(hst->GetMinimumBin());
      float mxh = hst->GetBinContent(hst->GetMaximumBin());
      if (mn > mnh)
        mn = mnh;
      if (mx < mxh)
        mx = mxh;
      if (!hbase)
        hbase = hst;
    }
  }
  return hbase;
}

void HistoManager::setHistosMinMaxRange(TVirtualPad* pad, float mn, float mx, float marginH, float marginL)
{
  float delta = mx - mn;
  if (delta <= 0)
    return;
  if (marginH > 0)
    mx += delta * marginH;
  else
    mx -= marginH;
  if (marginL > 0)
    mn -= delta * marginL;
  else
    mn += marginL;

  if (!pad)
    pad = (TPad*)gPad;
  if (!pad)
    return;
  TH1* hbase = 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TH1* hst = 0;
  for (int i = 0; i < size; i++) {
    TObject* obj = lst->At(i);
    if (!obj)
      continue;
    if (obj->InheritsFrom("TH1")) {
      hst = (TH1*)obj;
      hst->SetMinimum(mn);
      hst->SetMaximum(mx);
    }
  }
  pad->Modified();
  pad->Update();
}

TFrame* HistoManager::getFrame(TPad* pad)
{
  if (!pad)
    pad = (TPad*)gPad;
  if (!pad)
    return 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TFrame* frm = 0;
  for (int i = 0; i < size; i++) {
    TObject* obj = lst->At(i);
    if (!obj)
      continue;
    if (obj->InheritsFrom("TFrame")) {
      frm = (TFrame*)obj;
      break;
    }
  }
  return frm;
}

TPaveStats* HistoManager::getStatPad(TH1* hst)
{
  TList* lst = hst->GetListOfFunctions();
  if (!lst)
    return 0;
  int nf = lst->GetSize();
  for (int i = 0; i < nf; i++) {
    TPaveStats* fnc = (TPaveStats*)lst->At(i);
    if (fnc->InheritsFrom("TPaveStats"))
      return fnc;
  }
  return 0;
  //
}

TPaveStats* HistoManager::setStatPad(TH1* hst, float x1, float x2, float y1, float y2, Int_t stl, Int_t col)
{
  TPaveStats* pad = getStatPad(hst);
  if (!pad)
    return 0;
  pad->SetX1NDC(x1);
  pad->SetX2NDC(x2);
  pad->SetY1NDC(y1);
  pad->SetY2NDC(y2);
  if (stl >= 0)
    pad->SetFillStyle(stl);
  if (col >= 0)
    pad->SetTextColor(col);
  //
  gPad->Modified();
  return pad;
}

void HistoManager::setHStyle(TH1* hst, int col, int mark, float mrsize)
{
  hst->SetLineColor(col);
  hst->SetMarkerColor(col);
  //  hst->SetFillColor(col);
  hst->SetMarkerStyle(mark);
  hst->SetMarkerSize(mrsize);
  TList* lst = hst->GetListOfFunctions();
  if (lst) {
    int nf = lst->GetSize();
    for (int i = 0; i < nf; i++) {
      TObject* fnc = lst->At(i);
      if (fnc->InheritsFrom("TF1")) {
        ((TF1*)fnc)->SetLineColor(col);
        ((TF1*)fnc)->SetLineWidth(1);
        ((TF1*)fnc)->ResetBit(TF1::kNotDraw);
      } else if (fnc->InheritsFrom("TPaveStats")) {
        ((TPaveStats*)fnc)->SetTextColor(col);
      }
    }
  }
}

void HistoManager::setGStyle(TGraph* hst, int col, int mark, float mrsize)
{
  hst->SetLineColor(col);
  hst->SetMarkerColor(col);
  hst->SetFillColor(col);
  hst->SetMarkerStyle(mark);
  hst->SetMarkerSize(mrsize);
  TList* lst = hst->GetListOfFunctions();
  if (lst) {
    int nf = lst->GetSize();
    for (int i = 0; i < nf; i++) {
      TObject* fnc = lst->At(i);
      if (fnc->InheritsFrom("TF1")) {
        ((TF1*)fnc)->SetLineColor(col);
        ((TF1*)fnc)->SetLineWidth(1);
        ((TF1*)fnc)->ResetBit(TF1::kNotDraw);
      } else if (fnc->InheritsFrom("TPaveStats")) {
        ((TPaveStats*)fnc)->SetTextColor(col);
      }
    }
  }
}

//_________________________________________________________________________
TH1* HistoManager::prof2TH1(TProfile* prf, const char* addName)
{
  // convert profile to TH1 histo
  TString nm = prf->GetName();
  nm += addName;
  TAxis* xax = prf->GetXaxis();
  TH1* prof = 0;
  const TArrayD* bins = xax->GetXbins();
  if (bins->GetSize() == 0) {
    prof = new TH1F(nm.Data(), nm.Data(), prf->GetNbinsX(), xax->GetXmin(), xax->GetXmax());
  } else {
    prof = new TH1F(nm.Data(), nm.Data(), xax->GetNbins(), bins->GetArray());
  }
  for (int i = 1; i <= prof->GetNbinsX(); i++) {
    prof->SetBinContent(i, prf->GetBinContent(i));
    prof->SetBinError(i, prf->GetBinError(i));
  }
  return prof;
}

//_________________________________________________________________________
TH1* HistoManager::invertHisto(TH1* h)
{
  int nb = h->GetNbinsX();
  TH1* hi = (TH1*)h->Clone(Form("%s_inv", h->GetName()));
  for (int i = 0; i <= nb + 1; i++) {
    double v = h->GetBinContent(i);
    double e = h->GetBinError(i);
    int j = nb - i + 1;
    hi->SetBinContent(j, v);
    hi->GetBinError(j, e);
  }
  hi->SetTitle(Form("%s inv", h->GetTitle()));
  return hi;
}

TLatex* HistoManager::addLabel(const char* txt, float x, float y, int color, float size)
{
  TLatex* lt = new TLatex(x, y, txt);
  lt->SetNDC();
  lt->SetTextColor(color);
  lt->SetTextSize(size);
  lt->Draw();
  return lt;
}

TLegendEntry* addLegendEntry(TLegend* leg, const TH1* obj, const char* label, Option_t* option)
{
  if (!leg) {
    return nullptr;
  }
  auto ent = leg->AddEntry(obj, label, option);
  ent->SetTextColor(obj->GetLineColor());
  return ent;
}

TLegendEntry* addLegendEntry(TLegend* leg, const TGraph* obj, const char* label, Option_t* option)
{
  if (!leg) {
    return nullptr;
  }
  auto ent = leg->AddEntry(obj, label, option);
  ent->SetTextColor(obj->GetLineColor());
  return ent;
}

} // namespace na6p
