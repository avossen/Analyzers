import NovelFitters.NovelBaseFitter;
import NovelFitters.LundPID;
import NovelFitters.MyParticle;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;

import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;

import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.*;
import static java.nio.file.FileVisitResult.*;
import static java.nio.file.FileVisitOption.*;
import java.util.*;
import java.util.Arrays;
//import org.jlab.clas12.physics.*;

public class SimpleAnalyzer {

	protected H1F hLambdaMass;
	protected H1F hLambdaMassXfCut;
	protected H1F hLambdaXf;
	protected H1F hTrueXf;
	protected H2F hXfVsZ;
	
	protected H2F hPid2[];
	protected H1F hPid1[];

	protected H1F hhpid1;
	protected H1F hLambdaMassRes;
	protected H1F hTrueLambdaMass;
	protected H1F hLambdaMassXfCutRes;
	protected H1F hLambdaXfRes;
	protected H2F hXfVsZRes;
	public boolean isMC;

	protected H1F hPiPiMass;
	protected H1F hPiPiMassXfCut;
	protected H1F hPiPiXf;
	protected H2F hPiPiXfVsZ;

	protected H1F hDiPionMass;
	protected H1F hDiPionTheta;
	protected H1F hDiPionPPerp;
	protected H1F hDiPionPPerpLab2;
	
	protected H1F hDiPionMass2;
	protected H1F hDiPionTheta2;
	protected H1F hDiPionPPerp2;
	protected H1F hDiPionPhiH;
	protected H1F hDiPionPhiR;
	protected H1F hDiPionMissingMass;
	
	protected H1F hX;
	protected H1F hQ2;
	protected H1F hW;
	protected H1F hZ;
	protected H2F hQvsX;
	
	protected HipoDataSource reader;
	protected int matchCounter;
	protected int noMatchCounter;
	protected int m_numGoodFilterEvts;	
	protected int m_numEvtsWithPIDChi2;
	protected int m_numEvtsWithKinCuts;
	protected int m_numEventsWithBanks;
	//if one only wants to know how many events have at least one pair in it
	protected boolean m_EvtCountedKinCuts;
	protected boolean m_EvtCountedPIDChi2;
	protected boolean evtFulfillsMissingMass;

	protected final double m_pi = 0.1396;
	protected NovelBaseFitter novel_fitter;
	protected NovelBaseFitter novel_fitterMC;
	protected EventFilter filter;
	protected boolean printDebug=false;
	AsymData m_asymData;
	//ArrayList<EventData> m_eventData;
	//ArrayList<EventData> m_mcEventData;
	//hold reference to the current event
	EventData currentEvent;
	EventData currentMCEvent;

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			// exits program if input directory not specified
			System.out.println("ERROR: Please enter a hipo file as the first argument");
			System.exit(0);
		}
		SimpleAnalyzer analyzer = new SimpleAnalyzer();
		//analyzer.isMC=true;
		analyzer.isMC=false;
		analyzer.analyze(args);
		analyzer.plot();
	}
	
	//add the hadron pair to the event data
	protected HadronPairData addHadronPair(HadronPair pair, boolean isMC)
	{
		HadronPairData data=new HadronPairData();
		data.phiH=(float)pair.getPhiH();
		data.M=(float)pair.getMass();
		data.phiR=(float)pair.getPhiR();
		data.pTBreit=(float)pair.getPt();
		data.pTLab=(float)pair.getPtLab();
		data.xF=(float)pair.getXf();
		data.z=(float)pair.getZ();
		data.hasMC=false;
		data.theta=(float)pair.getTheta();
		data.theta1=(float)pair.getTheta1();
		data.theta2=(float)pair.getTheta2();
		data.phi1=(float)pair.getPhi1();
		data.phi2=(float)pair.getPhi2();
		data.mom1=(float)pair.getMom1();
		data.mom2=(float)pair.getMom1();
		
		
	//	System.out.println("saving, isMC?: " + isMC);
		
		//not the MC, so lets see if it has a matching MC that we want to save
		//the assumption is that java checks the references, so that things are not saved twice, just the reference
		if(!isMC)
		{
			if(pair.hasMatchingMC)
			{
				//System.out.println("has mc");
				data.hasMC=true;
				data.matchingMCPair=pair.matchingMCData;
			}
			this.currentEvent.pairData.add(data);
		}
		else
		{
			this.currentMCEvent.pairData.add(data);
		}
		return data;
	}
	
	public void analyze(String[] args) {
		reader= new HipoDataSource();
		//for debugging, use eventbuilder
		NovelBaseFitter.useStefanElectronCuts=false;
		NovelBaseFitter.useStefanHadronCuts=false;
		//for haruts mc
		NovelBaseFitter.useTimeBasedTracks=true;
		
		
		m_numGoodFilterEvts=0;
		m_numEvtsWithPIDChi2=0;
		m_numEvtsWithKinCuts=0;
		m_numEventsWithBanks=0;
		hhpid1=new H1F("hhpid1","hhpid1",30,-20,20);
		hPid1=new H1F[7];
		hPid2=new H2F[7];
		for(int i=0;i<7;i++)
		{
			String s="hpid"+i;
			hPid2[i]=new H2F(s,s,30,0,10,30,-20,20);
			s="hpid1"+i;
			hPid1[i]=new H1F(s,s,30,-20.0,20.0);
		}
		
		hLambdaMass = new H1F("lambdaMass", "lambdaMass", 100, 1.0, 1.2);
		hTrueLambdaMass = new H1F("trueLambdaMass", "trueLambdaMass", 100, 1.0, 2.0);
		hLambdaMassXfCut = new H1F("lambdaMass", "lambdaMass", 100, 1.0, 2.0);
		hLambdaXf = new H1F("Xf", "Xf", 100, -1.0, 1.0);
		hTrueXf = new H1F("trueXf", "trueXf", 100, -1.0, 1.0);
		hXfVsZ = new H2F("xfVsZ", "xfVsZ", 20, 0.0, 1.0, 20, 0.0, 1.0);
		
		hZ=new H1F("z","z",100,0.0,1.0);
	
		hQ2 = new H1F("Q2", "Q2", 100, 0.0, 10.0);
		hX = new H1F("x", "x", 100, 0.0, 1.0);
		hW = new H1F("W", "W", 100, 0.0, 5.0);
		hQvsX = new H2F("qVsX", "qVsX", 20, 0.0, 1.0, 20, 0.0, 10.0);
		

		hLambdaMassRes = new H1F("lambdaMassRes", "lambdaMassRes", 100, -0.5, 0.5);
		hLambdaMassXfCutRes = new H1F("lambdaMassRes", "lambdaMassRes", 100, -0.5, 0.5);
		hLambdaXfRes = new H1F("XfRes", "XfRes", 100, -0.3, 0.3);
		hXfVsZRes = new H2F("xfVsZResf", "xfVsZRes", 20, 0.0, 1.0, 20, 0.0, 1.0);

		hPiPiMass = new H1F("piPiMass", "piPiMass", 100, 0.3, 2.0);
		hPiPiMassXfCut = new H1F("piPiMass", "piPiMass", 100, 0.3, 2.0);
		hPiPiXf = new H1F("PiPiXf", "PiPiXf", 100, -1.0, 1.0);
		hPiPiXfVsZ = new H2F("piPixfVsZ", "PiPixfVsZ", 20, 0.0, 1.0, 20, 0.0, 1.0);

		hDiPionMass = new H1F("diPionMass", "diPionMass", 100, 0.0, 2.0);
		hDiPionTheta = new H1F("diPionTheta", "diPionTheta", 100, (-1) * Math.PI, Math.PI);
		hDiPionPPerp = new H1F("diPionPPerp", "diPionPPerp", 100, 0.0, 2.0);
		hDiPionPPerpLab2 = new H1F("diPionPPerpLab2", "diPionPPerpLab2", 100, 0.0, 3.0);
		
		hDiPionMass2 = new H1F("diPionMass2", "diPionMass2", 100, 0.0, 2.0);
		hDiPionTheta2 = new H1F("diPionTheta2", "diPionTheta2", 100, 0.0, Math.PI);
		hDiPionPPerp2 = new H1F("diPionPPerp2", "diPionPPerp2", 100, 0.0, 2.0);
		
		hDiPionPhiH=new H1F("diPIonPhiH","diPionPhiH",100,0.0,2*Math.PI);
		hDiPionPhiR=new H1F("diPIonPhiR","diPionPhiR",100,0.0,2*Math.PI);
		
		hDiPionMissingMass= new H1F("diPIonMissingMass","diPionMissingMass",100,0,6.0);

		int numLambdas=0;
		
		// define fitter class, argument is the beam energy (GeV)
		 novel_fitter = new NovelBaseFitter(10.6,false,false);
		//novel_fitter = new NovelBaseFitter(10.6, true, true);
		 if(isMC)
			 novel_fitterMC = new NovelBaseFitter(10.6, true, true);
		// define filter to apply to your events
		// here we look for events with one electron (11), one photon (22) (change to no
		// photon) and any number of other
		// positively charged particles (X+), negatively charged particles (X-) or
		// neutral
		// particles (Xn)
		//EventFilter filter = new EventFilter("11:X+:X-:Xn");
		 filter = new EventFilter("11:+211:-211:X+:X-:Xn");
		//check if args[0] is already a hipo file, else run over all in the folder
		Path singlefilename = Paths.get(args[0]);
		//two * crosses directory boundaries...
		PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:**.{hipo}");
		if (matcher.matches(singlefilename)) 
		{
			//System.out.println("matched single file");
			//singlefilename should already contain the full path, so don't add args[0], instead args[0] is the filename
			readData("", args[0]);
		}
		else
		{
			File folder = new File(args[0]);
			File[] listOfFiles = folder.listFiles();
			for (int iF = 0; iF < listOfFiles.length; iF++) {
				if (listOfFiles[iF].isFile()) {
					System.out.println("File " + listOfFiles[iF].getName());
					Path filename = Paths.get(listOfFiles[iF].getName());
					if (matcher.matches(filename)) {
						readData(args[0], listOfFiles[iF].getName());
					}
				}
			}
		}
		System.out.println("num evts: " + this.m_numEventsWithBanks + " filtered: " + this.m_numGoodFilterEvts + " numWithQ2, W cuts: "+ this.m_numEvtsWithKinCuts+ " and with chi2pid cuts: "+ this.m_numEvtsWithPIDChi2);
		System.out.println("found " + numLambdas + " Lambdas");
	//	System.out.println("match found " + matchCounter + " nomatch: " + noMatchCounter + " percent "+ matchCounter/(matchCounter+noMatchCounter));
	}

	void readData(String args0, String filename)
	{
			//one asymmetry file per data input file
			this.m_asymData=new AsymData();
			System.out.println("matched" + filename);
			reader.open(args0 + filename); // open hipo file
			while (reader.hasEvent() == true) { // cycle through events
				// load next event in the hipo file
				//System.out.println("new event----\n\n");
				HipoDataEvent event = (HipoDataEvent) reader.getNextEvent();
				this.m_EvtCountedPIDChi2=false;
				this.m_EvtCountedKinCuts=false;
				// apply fitter to your data event
				
				PhysicsEvent generic_Event = novel_fitter.getPhysicsEvent(event);
				PhysicsEvent generic_EventMC=new PhysicsEvent();
				if(isMC)
				{
					generic_EventMC = novel_fitterMC.getPhysicsEvent(event);
				}
				// System.out.println("Q2: " + novel_fitter.getQ2() + " W: " +
				// novel_fitter.Walt);
					
				this.m_numEventsWithBanks++;
				
				//numLambdas+=novel_fitter.getNumLambda();
				//System.out.println("found2 " + novel_fitter.getNumLambda());
				if (filter.isValid(generic_Event) == true) { // apply filter to current event
					// look at all particles
					//System.out.println("valid event");
					if(printDebug)
					{
					   printEventInfo(generic_Event);
						
					}
					
					
					
					hQ2.fill(novel_fitter.getQ2());
					hX.fill(novel_fitter.getX());
					hW.fill(novel_fitter.getW());
					hQvsX.fill(novel_fitter.getX(), novel_fitter.getQ2());
					m_numGoodFilterEvts++;
					//System.out.println("Q2: " + novel_fitter.getQ2() + " w: " + novel_fitter.Walt);
					if (novel_fitter.getQ2() < 1.0 || novel_fitter.Walt < 2.0)
						continue;
					
					// grab electron
				//	System.out.println("q2 w good");
					Particle electron = generic_Event.getParticle("[11]");
				//	System.out.println("found electron with px: "+ electron.px());
					// get all Pions and then loop over them to get xF, Q2, theta etc
					if(isMC)
					{
						//-->
						associateMCWithData(generic_Event, generic_EventMC);
					}
					this.evtFulfillsMissingMass=false; 
					this.currentEvent=new EventData();
					currentEvent.y=(float)novel_fitter.getY();
					//System.out.println("y: " + novel_fitter.getY() + " q2 "+ novel_fitter.getQ2());
					currentEvent.Q2=(float)novel_fitter.getQ2();
					currentEvent.W=(float)novel_fitter.getW();
					currentEvent.x=(float)novel_fitter.getX();
					//need to get helicity, run and event number from MC event
					//because they are not stored in REC::Event record
					if(this.isMC)
					{
						currentEvent.beamHelicity=this.novel_fitterMC.getBeamHelicity();
					}
					else
					{
						currentEvent.beamHelicity=novel_fitter.getBeamHelicity();
					}
					this.currentMCEvent=new EventData();
					currentMCEvent.beamHelicity=0;
					doDiHadrons(generic_Event,generic_EventMC,novel_fitter,novel_fitterMC);
					//System.out.println("pair data size: "+currentEvent.pairData.size());
					if(printDebug)
					{
						printDiHadrons(generic_Event);
					}
					if(this.currentEvent.pairData.size()>0)
					{		
						//System.out.println("pair data size: "+currentEvent.pairData.size());
						m_asymData.eventData.add(this.currentEvent);
						m_asymData.eventDataMC.add(this.currentMCEvent) ;
					}
						
					//doLambdas(generic_Event,generic_EventMC,novel_fitter,novel_fitterMC);
					/**
					 * ArrayList<Particle> pions = (ArrayList<Particle>)
					 * generic_Event.getParticleListByPid(211); ArrayList<Particle> pionsMinus =
					 * (ArrayList<Particle>) generic_Event.getParticleListByPid(-211);
					 * System.out.println("found " + pions.size() + "pions+ and " +
					 * pionsMinus.size() +" pi-");
					 * 
					 * for(Particle pion: pions) {
					 * 
					 * //also try to find protons and combine to Lambdas List<Particle> protons =
					 * generic_Event.getParticlesByPid(2212); for(Particle proton: protons) {
					 * 
					 * 
					 * } }
					 */

					
				}
			}
			this.saveData(filename);
		    reader.close();
	}
	
	void doLambdas(PhysicsEvent generic_Event, PhysicsEvent generic_EventMC, NovelBaseFitter m_novel_fitter, NovelBaseFitter m_novel_fitterMC)
	{
		for (int i = 0; i < generic_Event.count(); i++) {
			
			MyParticle part = (MyParticle) generic_Event.getParticle(i);
			if(part.charge()!=0)
			{
				if(part.matchingMCPartIndex<0)
				{
					noMatchCounter++;
				}
				else
				{
					matchCounter++;
				}
			}
		
			// System.out.println("matching mc particle index: " +
			// part.matchingMCPartIndex);

			// System.out.println("found particle with pid: " +part.pid() + " energy: "+
			// part.e() + " theta: "+part.theta()/Math.PI *180 + "p: "+part.px() + ", " +
			// part.py() + ", " +part.pz());
			if (part.pid() == 2212) {

				// System.out.println("found proton");
				for (int j = 0; j < generic_Event.count(); j++) {
					MyParticle part2 = (MyParticle) generic_Event.getParticle(j);
					// Systefm.out.println("lookign at pid " + part2.pid());
					if (part2.pid() == (-211)) {
						// test K_S hypothesis by assuming that
						LorentzVector pionHyp = new LorentzVector();
						pionHyp.setPxPyPzM(part.px(), part.py(), part.pz(), m_pi);
						LorentzVector lambdaCandidate = new LorentzVector(part.px() + part2.px(),
								part.py() + part2.py(), part.pz() + part2.pz(),
								part.e() + part2.e());
						LorentzVector piPiCandidate = new LorentzVector(pionHyp.px() + part2.px(),
								pionHyp.py() + part2.py(), pionHyp.pz() + part2.pz(),
								pionHyp.e() + part2.e());

						hLambdaMass.fill(lambdaCandidate.mass());
						hPiPiMass.fill(piPiCandidate.mass());

						LorentzVector lvTarget = new LorentzVector();

						double z = lambdaCandidate.e() / novel_fitter.getq().e();

						double zPiPi = piPiCandidate.e() / novel_fitter.getq().e();
						LorentzVector boostedLambda = new LorentzVector(lambdaCandidate);
						LorentzVector boostedPiPi = new LorentzVector(piPiCandidate);
						boostedLambda.boost(novel_fitter.gNBoost);
						boostedPiPi.boost(novel_fitter.gNBoost);
						double xF = boostedLambda.pz() / novel_fitter.Walt;
						double xFPiPi = boostedLambda.pz() / novel_fitter.Walt;
						if (part.matchingMCPartIndex != -1 && part2.matchingMCPartIndex != -1) {
							System.out.println("found lambda candidate with matching MC!!");

							Particle mc1 = generic_EventMC.getParticle(part.matchingMCPartIndex);
							Particle mc2 = generic_EventMC.getParticle(part2.matchingMCPartIndex);
							LorentzVector mcTruth = new LorentzVector(mc1.px() + mc2.px(),
							mc1.py() + mc2.py(), mc1.pz() + mc2.pz(), mc1.e() + mc2.e());
							hLambdaMassRes.fill(mcTruth.mass() - lambdaCandidate.mass());
							LorentzVector boostedMCTruth = new LorentzVector(mcTruth);
							boostedMCTruth.boost(novel_fitter.gNBoost);
							double xFTrue = boostedMCTruth.pz() / novel_fitter.Walt;
							hLambdaXfRes.fill(xFTrue - xF);
							hTrueLambdaMass.fill(mcTruth.mass());
							if (mcTruth.mass() < 1.12 && mcTruth.mass() > 1.11)
								hTrueXf.fill(xFTrue);
							// if(xFTrue>0.0)
							{
								// hLambdaMassXfCutRes.fill(mcTruth.mass()-lambdaCandidate.mass());
							}
						}

						hLambdaXf.fill(xF);
						hPiPiXf.fill(xF);
						if (xF > 0.0) {
							hLambdaMassXfCut.fill(lambdaCandidate.mass());
							if (lambdaCandidate.mass() < 1.25)
								hXfVsZ.fill(z, xF);
							// System.out.println("walt: " + novel_fitter.Walt + "boosted pz
							// "+boostedLambda.pz());
							// System.out.println("found lambda candidate with mass:
							// "+lambdaCandidate.mass()+ " xF: "+ xF);
						}
						if (xFPiPi > 0.0) {
							hPiPiMassXfCut.fill(lambdaCandidate.mass());
							if (piPiCandidate.mass() < 1.25)
								hPiPiXfVsZ.fill(z, xF);
							// System.out.println("walt: " + novel_fitter.Walt + "boosted pz
							// "+boostedLambda.pz());
							// System.out.println("found lambda candidate with mass:
							// "+lambdaCandidate.mass()+ " xF: "+ xF);
						}	
					}
				}
			}

		}
	
	}
	
	void doDiHadrons(PhysicsEvent generic_Event, PhysicsEvent generic_EventMC, NovelBaseFitter m_novel_fitter, NovelBaseFitter m_novel_fitterMC) {
		//System.out.println("in do dihad with "+generic_Event.count() + "particles ");
		for (int i = 0; i < generic_Event.count(); i++) 
		
		
		{
			MyParticle part = (MyParticle) generic_Event.getParticle(i);
			int sec= part.FTOFsector;
			if(sec<=0)
				sec=6;
			
			//System.out.println("sec: " + sec+ ", beta is " + part.beta);
			if(part.beta<10.0 && part.beta>-10.0)
			{
				this.hPid1[sec].fill(part.beta);
				this.hPid2[sec].fill(part.p(),part.beta);
				
				
				
				
				if(sec!=6)
				{
					//System.out.println("fill with beta "+ part.beta);
				    this.hhpid1.fill((float) part.beta);
				}
			}
			//System.out.println("time is: " + part.FTOFTime + " sector: " + part.FTOFsector);
			// System.out.println("matching mc particle index: " +
			// part.matchingMCPartIndex);

			
			
			if (part.pid() == LundPID.Pion.lundCode() || part.pid()==LundPID.Kaon.lundCode()) {
				
				for (int j = 0; j < generic_Event.count(); j++) {
					MyParticle part2 = (MyParticle) generic_Event.getParticle(j);
					// Systefm.out.println("lookign at pid " + part2.pid());
					if (part2.pid() == ((-1)*LundPID.Pion.lundCode()) || part2.pid()==((-1)*LundPID.Kaon.lundCode())){
						
						if(!this.m_EvtCountedKinCuts)
						{
							this.m_numEvtsWithKinCuts++;
							this.m_EvtCountedKinCuts=true;
						}
					
						
						HadronPair pair=new HadronPair(part,part2,m_novel_fitter.getq(),m_novel_fitter.getL(),m_novel_fitter.Walt,m_novel_fitter.gNBoost);

							
						
						hDiPionMass2.fill(pair.getMass());
						 hDiPionTheta2.fill(pair.getTheta());
						 hDiPionPPerp2.fill(pair.getPt());
						 hDiPionPPerpLab2.fill(pair.getPtLab());
						 hDiPionPhiH.fill(pair.getPhiH());
						 hDiPionPhiR.fill(pair.getPhiR());
						 hZ.fill(pair.getZ());
						 hDiPionMissingMass.fill(pair.getMissingMass());
						 
						if(pair.getMissingMass()>1.05 && pair.getZ()<0.95)
						{
							evtFulfillsMissingMass=true;
							if(part.m_chi2pid<=5.0 && part2.m_chi2pid<=5.0)
							{
								if(false==this.m_EvtCountedPIDChi2)
								{
									this.m_numEvtsWithPIDChi2++;
									this.m_EvtCountedPIDChi2=true;
								}
							}
						}
						else
						{
							continue;
						}
					
						
						LorentzVector pionPair = new LorentzVector(part.px() + part2.px(), part.py() + part2.py(),
								part.pz() + part2.pz(), part.e() + part2.e());
						hDiPionMass.fill(pionPair.mass());
						LorentzVector boostedPair = new LorentzVector(pionPair);
						boostedPair.boost(m_novel_fitter.gNBoost);
						// we are in the breit frame now, so I guess that the transverse component
						// should already be Pht
						double px = boostedPair.vect().x();
						double py = boostedPair.vect().y();
						double phT = Math.sqrt(px * px + py * py);
						hDiPionPPerp.fill(phT);
						//check if there is a MC counterpart and save
						
						if (part.matchingMCPartIndex != -1 && part2.matchingMCPartIndex != -1) {
							//System.out.println("found di hadron  candidate with matching MC!!");

							MyParticle mc1 = (MyParticle)generic_EventMC.getParticle(part.matchingMCPartIndex);
							MyParticle mc2 = (MyParticle)generic_EventMC.getParticle(part2.matchingMCPartIndex);
							HadronPair pairMC=new HadronPair(mc1,mc2,m_novel_fitterMC.getq(),m_novel_fitterMC.getL(),m_novel_fitterMC.Walt,m_novel_fitterMC.gNBoost);
							HadronPairData hpd=addHadronPair(pairMC,true);
							//System.out.println("found matching hadron pair");
							pair.hasMatchingMC=true;
							pair.matchingMC=pairMC;
							pair.matchingMCData=hpd;
						}
					//System.out.println("adding hadron pair");	
					//System.out.println("pair data before: " + this.currentEvent.pairData.size());	
						addHadronPair(pair,false);	
					//	System.out.println("pair data now: " + this.currentEvent.pairData.size());	
						
					}
				}
			}
		}

	}

	
	public void saveData(String hipoFilename)
	{
		String filename=hipoFilename.substring(0, hipoFilename.lastIndexOf('.'));
		filename=filename.substring(filename.lastIndexOf('/')+1,filename.length());
		filename=filename+".srn";
		System.out.println("saving java output to: " + filename);
		try
		{
			FileOutputStream file=new FileOutputStream(filename);
			ObjectOutputStream out = new ObjectOutputStream(file);
			out.writeObject(this.m_asymData);
			out.close();
			file.close();
			System.out.println("saved java output to: " + filename);
		}
		catch(IOException ex)
		{
			System.out.println("IOException is caught");
		}
	}
	
	
	public void plot() {
		System.out.println("in plotting function");
		EmbeddedCanvas can_lambda = new EmbeddedCanvas();
		can_lambda.setSize(1200, 600);
		can_lambda.divide(2, 2);
		can_lambda.setAxisTitleSize(24);
		can_lambda.setAxisFontSize(24);
		can_lambda.setTitleSize(24);
		can_lambda.cd(0);
		can_lambda.draw(hLambdaMass);
		can_lambda.cd(1);
		can_lambda.draw(hLambdaXf);
		can_lambda.cd(2);
		can_lambda.draw(hLambdaMassXfCut);
		can_lambda.cd(3);
		can_lambda.draw(hXfVsZ);
		// can_lambda.cd(1);can_e_ecal.draw(H_ESampl_ECal);
		// can_lambda.draw(g_m_ESampl_ECal,"same");
		// can_lambda.draw(g_s_ESampl_ECal,"same");
		can_lambda.save("lambda.png");

		EmbeddedCanvas can_piPi = new EmbeddedCanvas();
		can_piPi.setSize(1200, 600);
		can_piPi.divide(2, 2);
		can_piPi.setAxisTitleSize(24);
		can_piPi.setAxisFontSize(24);
		can_piPi.setTitleSize(24);
		can_piPi.cd(0);
		can_piPi.draw(hPiPiMass);
		can_piPi.cd(1);
		can_piPi.draw(hPiPiXf);
		can_piPi.cd(2);
		can_piPi.draw(hPiPiMassXfCut);
		can_piPi.cd(3);
		can_piPi.draw(hPiPiXfVsZ);
		// can_piPi.cd(1);can_e_ecal.draw(H_ESampl_ECal);
		// can_piPi.draw(g_m_ESampl_ECal,"same");
		// can_piPi.draw(g_s_ESampl_ECal,"same");
		can_piPi.save("piPiCandidates.png");

		EmbeddedCanvas can_Res = new EmbeddedCanvas();
		can_Res.setSize(1200, 600);
		can_Res.divide(2, 2);
		can_Res.setAxisTitleSize(24);
		can_Res.setAxisFontSize(24);
		can_Res.setTitleSize(24);
		can_Res.cd(0);
		can_Res.draw(hLambdaMassRes);
		can_Res.cd(1);
		can_Res.draw(hLambdaXfRes);
		can_Res.cd(2);
		can_Res.draw(hTrueLambdaMass);
		can_Res.cd(3);
		can_Res.draw(hTrueXf);
		// can_piPi.cd(1);can_e_ecal.draw(H_ESampl_ECal);
		// can_piPi.draw(g_m_ESampl_ECal,"same");
		// can_piPi.draw(g_s_ESampl_ECal,"same");
		can_Res.save("ResolutionStudies.png");

		EmbeddedCanvas can_dihad = new EmbeddedCanvas();
		can_dihad.setSize(1200, 600);
		can_dihad.divide(2, 2);
		can_dihad.setAxisTitleSize(24);
		can_dihad.setAxisFontSize(24);
		can_dihad.setTitleSize(24);
		can_dihad.cd(0);
		can_dihad.draw(hDiPionMass);
		can_dihad.cd(1);
		can_dihad.draw(hDiPionPPerp);
		can_dihad.save("dihadrons.png");
		
		EmbeddedCanvas can_dihad2 = new EmbeddedCanvas();
		can_dihad2.setSize(1800, 1000);
		can_dihad2.divide(3, 3);
		can_dihad2.setAxisTitleSize(30);
		can_dihad2.setAxisFontSize(30);
		can_dihad2.setTitleSize(30);
		can_dihad2.cd(0);
		can_dihad2.draw(hDiPionMass2);
		can_dihad2.cd(1);
		can_dihad2.draw(hDiPionPPerp2);
		can_dihad2.cd(2);
		can_dihad2.draw(hDiPionPPerpLab2);
		
		can_dihad2.cd(3);
		can_dihad2.draw(hDiPionTheta2);
		can_dihad2.cd(4);
		can_dihad2.draw(hDiPionPhiH);
		can_dihad2.cd(5);
		can_dihad2.draw(hDiPionPhiR);
		can_dihad2.cd(6);
		can_dihad2.draw(this.hDiPionMissingMass);
		can_dihad2.save("dihadrons2.png");
		
		EmbeddedCanvas can_kinematics = new EmbeddedCanvas();
		can_kinematics.setSize(1200, 600);
		can_kinematics.divide(2, 2);
		can_kinematics.setAxisTitleSize(24);
		can_kinematics.setAxisFontSize(24);
		can_kinematics.setTitleSize(24);
		//can_kinematics.setl
		can_kinematics.cd(0);
		can_kinematics.draw(hQ2);
		can_kinematics.getPad(0).getAxisY().setLog(true);
		//hQ2.getAxis().
		can_kinematics.cd(1);
		can_kinematics.getPad(1).getAxisX().setLog(true);
		can_kinematics.draw(hX);
		can_kinematics.cd(2);
		can_kinematics.draw(hW);
		can_kinematics.cd(3);
		can_kinematics.getPad(3).getAxisX().setLog(true);
		//can_kinematics.getPad(3).getAxisY().setLog(true);
		can_kinematics.getPad(3).getAxisZ().setLog(true);
		can_kinematics.draw(hQvsX);
		can_kinematics.save("kinematics.png");
		
		EmbeddedCanvas can_pid1=new EmbeddedCanvas();
		can_pid1.setSize(1200,600);
		can_pid1.divide(4, 2);
		for(int i=0;i<7;i++)
		{
			can_pid1.cd(i);
			can_pid1.draw(this.hPid1[i]);	
		}
		can_pid1.cd(7);
		can_pid1.draw(hhpid1);
		can_pid1.save("pid1.png");
		System.out.println("saved pid1");
		
		EmbeddedCanvas can_pid2=new EmbeddedCanvas();
		can_pid2.setSize(1200,600);
		can_pid2.divide(4, 2);
		for(int i=0;i<7;i++)
		{
			can_pid2.cd(i);
			can_pid2.getPad(i).getAxisZ().setLog(true);
			
		//	System.out.println("aobut to draw pid2 ");
			can_pid2.draw(this.hPid2[i]);	
		}
		can_pid2.save("pid2.png");

	}

	protected void associateMCWithData(PhysicsEvent generic_Event, PhysicsEvent generic_EventMC) {
//System.out.println("associate");
		for (int i = 0; i < generic_Event.count(); i++) {
			double minMomDiff = 1000.0;
			int minMomDiffIndex = -1;
			MyParticle part = (MyParticle) generic_Event.getParticle(i);
			double px = part.px();
			double py = part.py();
			double pz = part.pz();
			// System.out.println("data part px " + px + " py " + py + " pz " + pz);

			double mom = Math.sqrt(px * px + py * py + pz * pz);
			double theta = Math.toDegrees(Math.atan2(Math.sqrt(px * px + py * py), pz));
			double phi = Math.toDegrees(Math.atan2(py, px));
			// System.out.println("data theta "+ theta + " phi " + phi + " mom " +mom);
			for (int j = 0; j < generic_EventMC.count(); j++) {
				MyParticle partMC = (MyParticle) generic_EventMC.getParticle(j);
				double pxMC = partMC.px();
				double pyMC = partMC.py();
				double pzMC = partMC.pz();

				double momMC = Math.sqrt(pxMC * pxMC + pyMC * pyMC + pzMC * pzMC);
				double thetaMC = Math.toDegrees(Math.atan2(Math.sqrt(pxMC * pxMC + pyMC * pyMC), pzMC));
				double phiMC = Math.toDegrees(Math.atan2(pyMC, pxMC));
				if (i == 0) {
					// System.out.println("MC part px " + pxMC + " py " + pyMC + " pz " + pzMC);
					// System.out.println("mc teta "+ thetaMC + " phi " + phiMC + " mom " +momMC);
				}

				// System.out.println("Looking at MC part "+ j + " rel momDiff: " + ((mom-momMC)/momMC));
			//	if (Math.abs(mom - momMC) < 0.025 * momMC && Math.abs(theta - thetaMC) < 1.0 && Math.abs(phi - phiMC) < 5)
				//let's increase this
				if (Math.abs(mom - momMC) < 0.05 * momMC && Math.abs(theta - thetaMC) < 2.0 && Math.abs(phi - phiMC) < 10)
				{
					if (Math.abs(mom - momMC) < minMomDiff) {
						minMomDiff = Math.abs(mom - momMC);
						minMomDiffIndex = j;
					}
				}
			}
		//	System.out.println("associate mc part "+minMomDiffIndex + " with "+ i);
			part.matchingMCPartIndex = minMomDiffIndex;
		}

	}
	void printEventInfo(PhysicsEvent generic_Event)
	{
		double Q2=novel_fitter.getQ2();
		double W=novel_fitter.getW();
		double x=novel_fitter.getX();
		double y=novel_fitter.getY();
		
		System.out.printf("Run: %d Evt %d: Q2: %.2f x: %.2f W: %.2f y:: %.2f\n", novel_fitter.getRunNumber(),novel_fitter.getEvtNumber(),Q2,x,W,y);
		System.out.println("Found " + (generic_Event.count()) + " particles: ");
		for (int i = 0; i < generic_Event.count(); i++) 
		{
			MyParticle part = (MyParticle) generic_Event.getParticle(i);
			int pid=part.pid();
			double mom=part.p();
			int charge=part.charge();
			System.out.printf("PID: %d p: %.2f\n", pid,mom);
		
		}
		
	}
	void printDiHadrons(PhysicsEvent generic_Event)
	{
		for (int i = 0; i < generic_Event.count(); i++) 
		{
			MyParticle part = (MyParticle) generic_Event.getParticle(i);
			//System.out.println("lookint at first pid " + part.pid());
			if (part.pid() == LundPID.Pion.lundCode()) {
			//	System.out.println("lookint at first pid " + part.pid());
				for (int j = 0; j < generic_Event.count(); j++) {
					MyParticle part2 = (MyParticle) generic_Event.getParticle(j);
			//		System.out.println("lookign at second pid " + part2.pid());
					if (part2.pid() == ((-1)*LundPID.Pion.lundCode())){
						
						HadronPair pair=new HadronPair(part,part2,novel_fitter.getq(),novel_fitter.getL(),novel_fitter.Walt,novel_fitter.gNBoost);
						double pMass=pair.getMass();
						double pTheta=pair.getTheta();
						double xF=pair.getXf();
						double pPhiR=pair.getPhiR();
						double pPhiH=pair.getPhiH();	
						double pPt=pair.getPt();
						System.out.printf("pair mass: %.2f pT: %.2f xF: %.2f, theta: %.2f, phiR: %.2f, phiH: %.2f\n", pMass,pPt,xF,pTheta,pPhiR,pPhiH);
					}
				}
			}
		}
	}

}
