// NA6PCCopyright
#ifndef NA6P_MCAPP_H_
#define NA6P_MCAPP_H_

#include <TVirtualMCApplication.h>
#include <NA6PDetector.h>

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
  
private:
  NA6PDetector mDet;
};

#endif
