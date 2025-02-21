// NA6PCCopyright
#ifndef NA6P_MCAPP_H_
#define NA6P_MCAPP_H_

#include <TVirtualMCApplication.h>
#include <TParticle.h>

class NA6PMCStack;
class NA6PGenerator;
class TFile;
class TTree;

class NA6PMC : public TVirtualMCApplication
{
 public:
  NA6PMC(const char* name, const char* title);
  virtual ~NA6PMC();

  void ConstructGeometry() override;
  void ConstructOpGeometry() override;
  void InitGeometry() override;
  void AddParticles() override;
  void GeneratePrimaries() override;
  void BeginEvent() override;
  void FinishEvent() override;
  void BeginPrimary() override;
  void FinishPrimary() override;
  void Stepping() override;
  void PreTrack() override {}
  void PostTrack() override {}

  bool setupGenerator(const std::string& s);
  auto getGenerator() const { return mGenerator.get(); }

  void setRandomSeed(Long64_t r);
  auto getRandomSeed() const { return mRandomSeed; }

  void setVerbosity(int v) { mVerbosity = v; }
  auto getVerbosity() const { return mVerbosity; }

  void init();

  void selectTracksToSave();
  void createKineOutput(const std::string& outDir);
  void closeKineOutput();
  void writeKine();
  void forceCharmHadronicDecays();

 private:
  void clearHits();
  void addSpecialParticles();

  std::unique_ptr<NA6PMCStack> mStack{};
  std::unique_ptr<NA6PGenerator> mGenerator{};

  std::vector<TParticle> mMCTracks, *mMCTracksPtr = &mMCTracks;
  std::vector<int> mRemap; // tmp vector for selected tracks remapping
  std::vector<int> mSavID; // tmp vector for selected tracks original indices
  std::vector<std::pair<int, int>> mDtList;
  int mVerbosity = 0;
  ULong64_t mRandomSeed = 0;
  size_t mEvCount = 0;
  TFile* mKineFile = nullptr;
  TTree* mKineTree = nullptr;
};

#endif
