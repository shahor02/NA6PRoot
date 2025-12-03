// NA6PCCopyright

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NA6PBaseCluster + ;
#pragma link C++ class std::vector < NA6PBaseCluster> + ;
#pragma link C++ class NA6PTrack + ;
#pragma link C++ class std::vector < NA6PTrack> + ;

#pragma link C++ class NA6PReconstruction + ;
#pragma link C++ class NA6PVerTelReconstruction + ;
#pragma link C++ class ExtTrackPar + ;
#pragma link C++ class NA6PTrackerCA + ;
#pragma link C++ class NA6PFastTrackFitter + ;

#pragma link C++ class NA6PRecoParam + ;
#pragma link C++ class na6p::conf::ConfigurableParamHelper < NA6PRecoParam> + ;

#endif
