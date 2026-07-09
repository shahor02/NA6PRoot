#include <cstdint>
#include "NA6PEventReader.h"

void readEvent(
    const char* fileNameVerTel         = "VerticesVerTel.root",
    const char* fileNameTracksVerTel   = "TracksVerTel.root",
    const char* fileNameTracksMuonSpec = "TracksMuonSpec.root",
    const char* fileNameTracksMatching = "TracksMatching.root",
    const char* fileNameMC             = "MCKine.root",
    bool readMC                        = true
)
{
  NA6PEventReader reader(fileNameVerTel, fileNameTracksVerTel, fileNameTracksMuonSpec, fileNameTracksMatching, fileNameMC, readMC);

  for (std::int64_t iev = 0; iev < reader.entries(); ++iev) {
    reader.loadEvent(iev);

    const auto& vertices = reader.verticesVerTel();
    const auto& vtTracks = reader.tracksVerTel();
    const auto& msTracks = reader.tracksMuonSpec();
    const auto& matches  = reader.matches();
    const auto& mcParts  = reader.mcParticles();

    std::cout << "event " << iev << " nVtx=" << vertices.size() << " nVT=" << vtTracks.size() << " nMS=" << msTracks.size() << " nMatch=" << matches.size() << " nMC=" << mcParts.size() << std::endl;
  }
}
