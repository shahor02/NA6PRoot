// NA6PCCopyright
//
// based on
//  marian.ivanov@cern.ch
//
//  ------------------------------------------------------------------------------------------------
//  TreeStream
//  Standard stream (cout) like input for the tree
//  Run and see TreeStreamer::Test() - to see TreeStreamer functionality
//  ------------------------------------------------------------------------------------------------
//
//  -------------------------------------------------------------------------------------------------
//  TreeSRedirector
//  Redirect file to  different TreeStreams
//  Run and see   TreeSRedirector::Test() as an example of TreeSRedirector functionality
//
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
#ifndef NA6P_TREESTREAM_H
#define NA6P_TREESTREAM_H

#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TFile.h>
class TObjArray;
class TDataType;

class NA6PTreeDataElement : public TNamed
{
  friend class NA6PTreeStream;

 public:
  NA6PTreeDataElement(Char_t type);
  NA6PTreeDataElement(TDataType* type);
  NA6PTreeDataElement(TClass* cl);
  void SetPointer(void* pointer) { fPointer = pointer; }
  Char_t GetType() const { return fType; }

 protected:
  NA6PTreeDataElement(const NA6PTreeDataElement& tde);
  NA6PTreeDataElement& operator=(const NA6PTreeDataElement& tde);

  Char_t fType;      // type of data element
  TDataType* fDType; // data type pointer
  TClass* fClass;    // data type pointer
  void* fPointer;    // pointer to element

  ClassDef(NA6PTreeDataElement, 1);
};

class NA6PTreeStream : public TNamed
{
  friend class NA6PTreeSRedirector;

 public:
  NA6PTreeStream(const char* treename, TTree* externalTree = NULL);
  ~NA6PTreeStream();
  void Close();
  static void Test();
  Int_t CheckIn(Char_t type, void* pointer);
  // Int_t CheckIn(const char *type, void *pointer);
  Int_t CheckIn(TObject* o);
  void BuildTree();
  void Fill();
  Double_t GetSize() { return fTree->GetZipBytes(); }
  NA6PTreeStream& Endl();
  //
  NA6PTreeStream& operator<<(Bool_t& b)
  {
    CheckIn('B', &b);
    return *this;
  }
  NA6PTreeStream& operator<<(Char_t& c)
  {
    CheckIn('B', &c);
    return *this;
  }
  NA6PTreeStream& operator<<(UChar_t& c)
  {
    CheckIn('b', &c);
    return *this;
  }
  NA6PTreeStream& operator<<(Short_t& h)
  {
    CheckIn('S', &h);
    return *this;
  }
  NA6PTreeStream& operator<<(UShort_t& h)
  {
    CheckIn('s', &h);
    return *this;
  }
  NA6PTreeStream& operator<<(Int_t& i)
  {
    CheckIn('I', &i);
    return *this;
  }
  NA6PTreeStream& operator<<(UInt_t& i)
  {
    CheckIn('i', &i);
    return *this;
  }
  NA6PTreeStream& operator<<(Long_t& l)
  {
    CheckIn('L', &l);
    return *this;
  }
  NA6PTreeStream& operator<<(ULong_t& l)
  {
    CheckIn('l', &l);
    return *this;
  }
  NA6PTreeStream& operator<<(Long64_t& l)
  {
    CheckIn('L', &l);
    return *this;
  }
  NA6PTreeStream& operator<<(ULong64_t& l)
  {
    CheckIn('l', &l);
    return *this;
  }
  NA6PTreeStream& operator<<(Float_t& f)
  {
    CheckIn('F', &f);
    return *this;
  }
  NA6PTreeStream& operator<<(Double_t& d)
  {
    CheckIn('D', &d);
    return *this;
  }
  NA6PTreeStream& operator<<(TObject* o)
  {
    CheckIn(o);
    return *this;
  }
  NA6PTreeStream& operator<<(const Char_t* name);
  TTree* GetTree() const { return fTree; }

 protected:
  //

  NA6PTreeStream(const NA6PTreeStream& ts);
  NA6PTreeStream& operator=(const NA6PTreeStream& ts);

  TObjArray* fElements;   // array of elements
  TObjArray* fBranches;   // pointers to branches
  TTree* fTree;           // data storage
  Int_t fCurrentIndex;    // index of current element
  Int_t fId;              // identifier of layout
  TString fNextName;      // name for next entry
  Int_t fNextNameCounter; // next name counter
  Int_t fStatus;          // status of the layout

  ClassDef(NA6PTreeStream, 1);
};

class NA6PTreeSRedirector : public TObject
{
 public:
  NA6PTreeSRedirector(const char* fname = "", const char* option = "update");
  virtual ~NA6PTreeSRedirector();
  void Close();
  static void Test();
  static void Test2();
  static void UnitTestSparse(Double_t scale, Int_t testEntries);
  static void UnitTest(Int_t testEntries = 5000);
  void StoreObject(TObject* object);
  TFile* GetFile() { return fDirectory->GetFile(); }
  TDirectory* GetDirectory() { return fDirectory; }
  virtual NA6PTreeStream& operator<<(Int_t id);
  virtual NA6PTreeStream& operator<<(const char* name);
  void SetDirectory(TDirectory* sfile);
  void SetFile(TFile* sfile) { SetDirectory(sfile); }
  void SetExternalTree(const char* name, TTree* externalTree);
  static void SetDisabled(Bool_t b = kTRUE) { fgDisabled = b; }
  static Bool_t IsDisabled() { return fgDisabled; }
  static void FixLeafNameBug(TTree* tree);

 private:
  NA6PTreeSRedirector(const NA6PTreeSRedirector& tsr);
  NA6PTreeSRedirector& operator=(const NA6PTreeSRedirector& tsr);

  TDirectory* fDirectory;   // file
  Bool_t fDirectoryOwner;   // do we own the directory?
  TObjArray* fDataLayouts;  // array of data layouts
  static Bool_t fgDisabled; // disable - do not open any files

  ClassDef(NA6PTreeSRedirector, 1);
};

#endif
