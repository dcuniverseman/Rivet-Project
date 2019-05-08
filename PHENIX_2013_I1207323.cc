// -*- C++ -*-
#include "Rivet/HeavyIonAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet
{
/// **********************************************************************************
/// This is a Rivet Analysis for Medium modification factors at 200 GeV
///   It plots over events in the 0-40% centrality range
///   Triggers particles are defined as
///   gamma particles with 5 GeV < pT < 9 GeV and hadrons from .5 GeV < pT < 7 GeV
///   associated hadron yields are determined as a function of zT=phT/pγT


	class PHENIX_2013_I1207323 : public HeavyIonAnalysis
	{

	private:
		const int N_CENT_TYPES = 1; //Only 1 bin
		const int CENT_TYPE_EDGES[N_CENT_TYPES][2] = {{0.0,40}; // Centrality is 0-40%

		// histogram count;
		Scatter2DPtr _h1dPhi[N_CENT_TYPES]; // delta phi, split by centrality

		// Number of trigger particles for each centrality over
		//  all events. Needed in finalize() for normalization
		unsigned long long int nTrigger[N_CENT_TYPES];

	public:
		/// Constructor
		PHENIX_2013_I1207323():HeavyIonAnalysis("PHENIX_2013_I1207323") {}

		void init()
		{

			addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod");

			//**** Trigger particle set ****
			//  This is code to trigger off of a
			//   specific type or particle in a certain range. There will ALWAYS be an absolute eta
			//   cut on particles. It may differ from your trigger and associated
			//   particles.
			const int pidPI0 = 111; // The PID for Pi0
			const int pidGAMMA = 22; // The PID for Gamma
			//I believe I have a trigger hadron, but idk what it is.

			// I apply the following cuts for trigger particles:
			//  3.1 < |eta| < 3.9  for all particles
			//  and it is either
			//      A Pi0 with pT between .12 and .16 GeV *I think*
			//  or  A gamma with pT between 5 and 9 GeV
			Cuts cutTrigger = Cuts::abseta < 3.9  or abseta > 3.1 && ((Cuts::pid == pidPI0 && Cuts::pt > .12 * GeV && Cuts::pt < .16 * GeV) || (Cuts::pid == pidGAMMA && Cuts::pt > 5.0 * GeV && Cuts::pt < 9.0 * GeV))
			FinalState fs(cutTrigger);
			declare(fs, "partTrigger");

			//**** Associated particle set ****
			// Charged Particles with cuts
			// Max pT for an associated particle not explicitly stated in paper, have to come back and make more explicit values.
			Cuts cutAssoc = Cuts::abseta < 1.0 && Cuts::pt > 1.2 * GeV && Cuts::pt < 20.0 * GeV;
			ChargedFinalState cfs(cutAssoc);
			declare(cfs, "partAssoc");

			//**** Book Histograms ****
			//All are scatter plots

			_hist_c1 = bookScatter2D("d01-x01-y01");
			_hist2_c1 = bookScatter2D("d02-x01-y01");
			_hist3_c1 = bookScatter2D("d03-x01-y01");
			_hist4_c1 = bookScatter2D("d04-x01-y01");
			_hist5_c1 = bookScatter2D("d05-x01-y01");
			_hist6_c1 = bookScatter2D("d06-x01-y01");
			_hist7_c1 = bookScatter2D("d07-x01-y01");
			_hist8_c1 = bookScatter2D("d08-x01-y01");
			_hist9_c1 = bookScatter2D("d09-x01-y01");
			_hist10_c1 = bookScatter2D("d10-x01-y01");
			_hist11_c1 = bookScatter2D("d11-x01-y01");
			_hist12_c1 = bookScatter2D("d12-x01-y01");
			_hist13_c1 = bookScatter2D("d13-x01-y01");
			_hist14_c1 = bookScatter2D("d14-x01-y01");
			_hist15_c1 = bookScatter2D("d15-x01-y01");
			_hist16_c1 = bookScatter2D("d16-x01-y01");
			_hist17_c1 = bookScatter2D("d17-x01-y01");
			_hist18_c1 = bookScatter2D("d18-x01-y01");
			_hist19_c1 = bookScatter2D("d19-x01-y01");
			_hist20_c1 = bookScatter2D("d20-x01-y01");
			//**** Initialize counters ****
			// Set number of events in counters to zero
			// These are used later for normalizing histograms

			for (int i = 0; i < N_CENT_TYPES; ++i) nTrigger[i] = 0;
		}

		/// Per event calculations
		///  Do your per event calculations and
		///  fill your histograms here
		void analyze(const Event& event)
		{
			// Get the centrality for each event
			const double c = centrality(event, "IPMethod");

			// The first 50 events (number from addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 50, "IPMethod") )
			//  will give a centrality outside of 0-100, specifically -1.0
			if ((c < 0.) || (c > 100.)) vetoEvent;

			
	
			int centralityIndex = -1; // -1 - not in plot range
									
			// Find which centrality range of interest this falls into
			for (int i = 0; i < N_CENT_TYPES; ++i) {
				if (c > CENT_TYPE_EDGES[i][0] && c <= CENT_TYPE_EDGES[i][1]) {
					centralityIndex = i;
					break; // since theres no centrality overlap, break the loop when you find the correct index
				}
			}
			// If not a centrality of interest, stop processing the event
			if (centralityIndex == -1) vetoEvent;


			const FinalState & fs = apply<FinalState>(event,"partTrigger");
			const ChargedFinalState & cfs = apply<ChargedFinalState>(event,"partAssoc");

			// Get particles from projection objects.
			Particles tracksTrigger = fs.particlesByPt();
			Particles tracksAssoc = cfs.particlesByPt();

			// Increase appropriate counter
			nTrigger[centralityIndex] += tracksTrigger.size();

			// foreach goes through each Particle in tracksTrigger, and references it in partTrigger.
			// calculating delta-phi between trigger particles and charged hadrons
			double deltaPhi;
			foreach (const Particle& partTrigger, tracksTrigger) {
				// Loop over all associated particles
				foreach (const Particle& partAssoc, tracksAssoc) {
					// Only include associated particles with pT less than the trigger
					if (partAssoc.pt() < partTrigger.pt()) {
						deltaPhi = partAssoc.phi() - partTrigger.phi();
						//   so make sure deltaPhi falls in this range
						// M_PI is part of the <math.h> header
						while (deltaPhi < 0) deltaPhi +=  M_PI/2;

						// Fill plots
						//  In this case use a 1 as the second parameter, the weight.
						//  With event weighting factored in, that number would be variable
						//  and nTrigger would need to factor in weighting
						_h1dPhi[centralityIndex]->fill(deltaPhi,1);
					}
				}
			}
		}

		//  normalize my histogram by dividing by
		//  the number of triggers.
		void finalize()
		{
			for (int i = 0; i < N_CENT_TYPES; ++i){
				// calculate the inverse of the number of triggers
				double scaleFactor = 1./(double)nTrigger[i];

				// scale by that factor
				_h1dPhi[i]->scaleW(scaleFactor);
			}
		}
	};

	DECLARE_RIVET_PLUGIN(PHENIX_2013_I1207323);
}
