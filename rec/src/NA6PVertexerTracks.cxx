// NA6PCCopyright
#include <NA6PVertexerTracks.h>

ClassImp(NA6PVertexerTracks)

void NA6PVertexerTracks::createTracksPool(const std::vector<NA6PTrack>& tracks)
{
  mTracksPool.clear();
  auto ntGlo = tracks.size();
  mTracksPool.reserve(ntGlo);
  for (uint32_t i = 0; i < ntGlo; i++) {
    NA6PTrack trc = tracks[i];
    if (!trc.propagateToDCABeamAxis(mBeamX, mBeamY, mMaxDCA))
      continue;
    auto& tvf = mTracksPool.emplace_back(trc);
    if (!tvf.isValid()) {
      mTracksPool.pop_back(); // discard bad track
      continue;
    }
  }
}
