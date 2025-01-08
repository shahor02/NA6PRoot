// NA6PCCopyright

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class NA6PMCStack + ;
#pragma link C++ class NA6PGenerator + ;
#pragma link C++ class NA6PGenCutParam + ;
#pragma link C++ class na6p::conf::ConfigurableParamHelper < NA6PGenCutParam> + ;

#pragma link C++ class NA6PBaseHit + ;
#pragma link C++ class std::vector<NA6PBaseHit > + ;

#pragma link C++ class NA6PVerTelHit + ;
#pragma link C++ class std::vector<NA6PVerTelHit > + ;

#pragma link C++ class NA6PMuonSpecHit + ;
#pragma link C++ class std::vector<NA6PMuonSpecHit > + ;

#pragma link C++ class std::unordered_map<std::string, float > + ;

#pragma link C++ class NA6PMCEventHeader + ;
#pragma link C++ class std::vector<NA6PMCEventHeader > + ;

#pragma link C++ class NA6PMCGenHeader + ;
#pragma link C++ class std::vector<NA6PMCGenHeader > + ;

#pragma link C++ class NA6PGenBox + ;

#pragma link C++ class std::vector<TParticle > + ;



#endif
