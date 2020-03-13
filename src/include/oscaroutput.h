#ifndef INCLUDE_OSCAROUTPUT_H_
#define INCLUDE_OSCAROUTPUT_H_

#include "smash/file.h"
#include "smash/oscaroutput.h"
#include "smash/outputparameters.h"
#include "smash/particles.h"

#include "gen.h"
#include "params.h"

void write_oscar_output() {
  // initialize Oscar output
  const smash::bf::path OutputPath = params::sSpectraDir;
  std::unique_ptr<smash::OutputInterface> OscarOutput =
       create_oscar_output("Oscar2013", "Particles", OutputPath,
                           smash::OutputParameters());

  for(int iev=0; iev<params::NEVENTS; iev++){
    // Create Particles oject which contains the full particle list of each event
    std::unique_ptr<smash::Particles> particles = smash::make_unique<smash::Particles>();

    for(int inpart=0; inpart<gen::npart[iev]; inpart++) {
      smash::ParticleData* data = gen::pList[iev][inpart];
      particles->insert(*data);
    }
    // Write output
    OscarOutput->at_eventend(*particles, iev, 0.0, false);
  }
}

#endif // INCLUDE_OSCAROUTPUT_H_
