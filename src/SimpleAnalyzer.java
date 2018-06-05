import NovelFitters.NovelBaseFitter;
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

//import org.jlab.clas12.physics.*;

public class SimpleAnalyzer {

	protected H1F hLambdaMass;
	protected H1F hLambdaMassXfCut;
	protected H1F hLambdaXf;
	protected H1F hTrueXf;
	protected H2F hXfVsZ;
	
	protected H1F hLambdaMassRes;
	protected H1F hTrueLambdaMass;
	protected H1F hLambdaMassXfCutRes;
	protected H1F hLambdaXfRes;
	protected H2F hXfVsZRes;
	
	protected H1F hPiPiMass;
	protected H1F hPiPiMassXfCut;
	protected H1F hPiPiXf;
	protected H2F hPiPiXfVsZ;
	
	protected H1F hDiPionMass;
	protected H1F hDiPionTheta;
	protected H1F hDiPionPPerp;
	
	
	protected final double m_pi=0.1396;
	protected NovelBaseFitter novel_fitter;
	protected NovelBaseFitter novel_fitterMC;
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			// exits program if input directory not specified 
        	System.out.println("ERROR: Please enter a hipo file as the first argument");
       		System.exit(0);
    	}
		SimpleAnalyzer analyzer=new SimpleAnalyzer();
		analyzer.analyze(args);
		analyzer.plot();
	}
		
	public void analyze(String[] args)
	{
		hLambdaMass=new H1F("lambdaMass","lambdaMass",100,1.0,1.2);	
		hTrueLambdaMass=new H1F("trueLambdaMass","trueLambdaMass",100,1.0,2.0);
		hLambdaMassXfCut=new H1F("lambdaMass","lambdaMass",100,1.0,2.0);	
		hLambdaXf=new H1F("Xf","Xf",100,-1.0,1.0);	
		hTrueXf=new H1F("trueXf","trueXf",100,-1.0,1.0);	
		hXfVsZ=new H2F("xfVsZ","xfVsZ",20,0.0,1.0,20,0.0,1.0);
		
		hLambdaMassRes=new H1F("lambdaMassRes","lambdaMassRes",100,-0.5,0.5);	
		hLambdaMassXfCutRes=new H1F("lambdaMassRes","lambdaMassRes",100,-0.5,0.5);	
		hLambdaXfRes=new H1F("XfRes","XfRes",100,-0.3,0.3);	
		hXfVsZRes=new H2F("xfVsZResf","xfVsZRes",20,0.0,1.0,20,0.0,1.0);
		
		hPiPiMass=new H1F("piPiMass","piPiMass",100,0.3,2.0);	
		hPiPiMassXfCut=new H1F("piPiMass","piPiMass",100,0.3,2.0);	
		hPiPiXf=new H1F("PiPiXf","PiPiXf",100,-1.0,1.0);	
		hPiPiXfVsZ=new H2F("piPixfVsZ","PiPixfVsZ",20,0.0,1.0,20,0.0,1.0);
		
		hDiPionMass=new H1F("diPionMass","diPionMass",100,0.0, 2.0);
		hDiPionTheta=new H1F("diPionTheta","diPionTheta",100,(-1)*Math.PI, Math.PI);
		hDiPionPPerp=new H1F("diPionPPerp","diPionPPerp",100,0.0, 2.0);
		
		HipoDataSource reader = new HipoDataSource();
		// define fitter class, argument is the beam energy (GeV)
       // novel_fitter = new NovelBaseFitter(10.6,false,false);
		 novel_fitter = new NovelBaseFitter(10.6,true,true);
        novel_fitterMC = new NovelBaseFitter(10.6,true,true);
        // define filter to apply to your events
        // here we look for events with one electron (11), one photon (22) (change to no photon) and any number of other
        // positively charged particles (X+), negatively charged particles (X-) or neutral 
        // particles (Xn)
        EventFilter filter = new EventFilter("11:X+:X-:Xn");

		File folder = new File(args[0]);
		File[] listOfFiles = folder.listFiles();
		for (int iF = 0; iF < listOfFiles.length; iF++) {
			if (listOfFiles[iF].isFile()) {
				System.out.println("File " + listOfFiles[iF].getName());
		        	PathMatcher matcher =
		        	FileSystems.getDefault().getPathMatcher("glob:*.{hipo}");

		        	Path filename = Paths.get(listOfFiles[iF].getName());
		        	if (matcher.matches(filename)) {
		        	    System.out.println("matched" + filename);
		   
		        	    
		reader.open(args[0]+listOfFiles[iF].getName()); // open hipo file

		while(reader.hasEvent()==true){ // cycle through events
			// load next event in the hipo file
	//		System.out.println("new event----\n\n");
        	HipoDataEvent event = (HipoDataEvent) reader.getNextEvent();

            // apply fitter to your data event
            PhysicsEvent  generic_Event  =novel_fitter.getPhysicsEvent(event);
            PhysicsEvent  generic_EventMC  =novel_fitterMC.getPhysicsEvent(event);
      //      System.out.println("Q2: " + novel_fitter.getQ2() + " W: " + novel_fitter.Walt);
            if(novel_fitter.getQ2()<1.0 || novel_fitter.Walt<2.5)
            		continue;		
            if(filter.isValid(generic_Event)==true){ // apply filter to current event
//look at all particles
            	
                // grab electron
                Particle electron = generic_Event.getParticle("[11]");
               // System.out.println("found electron with px: "+ electron.px());
                //get all Pions and then loop over them to get xF, Q2, theta etc
                associateMCWithData(generic_Event, generic_EventMC);
               doDiHadrons(generic_Event); 
                
                
                
                	for(int i=0;i<generic_Event.count();i++)
                	{
                		MyParticle part=(MyParticle)generic_Event.getParticle(i);
    //            		System.out.println("matching mc particle index: " + part.matchingMCPartIndex);
                               		
                	//	System.out.println("found particle with pid: " +part.pid() + " energy: "+ part.e() + " theta: "+part.theta()/Math.PI *180 + "p: "+part.px() + ", " + part.py() + ", " +part.pz());     
                		if(part.pid()==2212)
                		{
                			
                		//	System.out.println("found proton");
                		 	for(int j=0;j<generic_Event.count();j++)
                        	{
                        		MyParticle part2=(MyParticle)generic_Event.getParticle(j);
                        	//	Systefm.out.println("lookign at pid " + part2.pid());
                        		if(part2.pid()==(-211))
                        		{
                        			//test K_S hypothesis by assuming that 
                        			LorentzVector pionHyp=new LorentzVector();
                        			pionHyp.setPxPyPzM(part.px(), part.py(), part.pz(), m_pi);
                        			LorentzVector lambdaCandidate=new LorentzVector(part.px()+part2.px(),part.py()+part2.py(),part.pz()+part2.pz(),part.e()+part2.e());
                        			LorentzVector piPiCandidate=new LorentzVector(pionHyp.px()+part2.px(),pionHyp.py()+part2.py(),pionHyp.pz()+part2.pz(),pionHyp.e()+part2.e());
                        	
                        			hLambdaMass.fill(lambdaCandidate.mass());
                        			hPiPiMass.fill(piPiCandidate.mass());
                        			
                        			
                        	
                        			
                        			LorentzVector lvTarget=new LorentzVector();
                        			
                        		
                        			double z=lambdaCandidate.e()/novel_fitter.getq().e();
                        			
                        			double zPiPi=piPiCandidate.e()/novel_fitter.getq().e();
                        			LorentzVector boostedLambda=new LorentzVector(lambdaCandidate);
                        			LorentzVector boostedPiPi=new LorentzVector(piPiCandidate);
                        			boostedLambda.boost(novel_fitter.gNBoost);
                        			boostedPiPi.boost(novel_fitter.gNBoost);
                        			double xF=boostedLambda.pz()/novel_fitter.Walt;
                        			double xFPiPi=boostedLambda.pz()/novel_fitter.Walt;
                        			if(part.matchingMCPartIndex!=-1 && part2.matchingMCPartIndex!=-1)
                        			{
                        				System.out.println("found lambda candidate with matching MC!!");
                        				
                        				Particle mc1=generic_EventMC.getParticle(part.matchingMCPartIndex);
                        				Particle mc2=generic_EventMC.getParticle(part2.matchingMCPartIndex);
                        				LorentzVector mcTruth=new LorentzVector(mc1.px()+mc2.px(),mc1.py()+mc2.py(),mc1.pz()+mc2.pz(),mc1.e()+mc2.e());
                        				hLambdaMassRes.fill(mcTruth.mass()-lambdaCandidate.mass());
                        				LorentzVector boostedMCTruth=new LorentzVector(mcTruth);
                        				boostedMCTruth.boost(novel_fitter.gNBoost);
                        				double xFTrue=boostedMCTruth.pz()/novel_fitter.Walt;
                        				hLambdaXfRes.fill(xFTrue-xF);
                        				hTrueLambdaMass.fill(mcTruth.mass());
                        				if(mcTruth.mass()<1.12 && mcTruth.mass()>1.11)
                        						hTrueXf.fill(xFTrue);
                        			//	if(xFTrue>0.0)
                        				{
                        				//	hLambdaMassXfCutRes.fill(mcTruth.mass()-lambdaCandidate.mass());	
                        				}
                        			}
                        				
                        			
                        			
                        			
                        			
                        			hLambdaXf.fill(xF);
                        			hPiPiXf.fill(xF);
                        			if(xF>0.0)
                        			{
                        				hLambdaMassXfCut.fill(lambdaCandidate.mass());
                        				if(lambdaCandidate.mass()<1.25)
                        					hXfVsZ.fill(z, xF);	
                        			//System.out.println("walt: " + novel_fitter.Walt + "boosted pz "+boostedLambda.pz());
                        			//System.out.println("found lambda candidate with mass: "+lambdaCandidate.mass()+ " xF: "+ xF);                      			                			
                        			}
                        			if(xFPiPi>0.0)
                        			{
                        				hPiPiMassXfCut.fill(lambdaCandidate.mass());
                        				if(piPiCandidate.mass()<1.25)
                        					hPiPiXfVsZ.fill(z, xF);	
                        			//System.out.println("walt: " + novel_fitter.Walt + "boosted pz "+boostedLambda.pz());
                        			//System.out.println("found lambda candidate with mass: "+lambdaCandidate.mass()+ " xF: "+ xF);                      			                			
                        			}
                        		}
                        	}
                		}
                		
                	}
                	
                /**
                ArrayList<Particle> pions = (ArrayList<Particle>) generic_Event.getParticleListByPid(211);
                ArrayList<Particle>	pionsMinus = (ArrayList<Particle>) generic_Event.getParticleListByPid(-211);
                System.out.println("found " + pions.size() + "pions+ and " + pionsMinus.size() +" pi-");
                
                for(Particle pion: pions)
                {
                	
                //also try to find protons and combine to Lambdas
                	   List<Particle> protons = generic_Event.getParticlesByPid(2212);
                       for(Particle proton: protons)
                       {
                       	
                       
                       }
                }
                */
                
                
                // grab first photon 
                //Particle photon = generic_Event.getParticle("[22]");

                // Particle class allows for addition of particles (adds four momenta)
                //Particle sum_particle = generic_Event.getParticle("[11]+[22]");
                // Particles have PID, vertices and momenta properties

                // print out the energy of the electron+photon
              //  System.out.println("The sum of the electron and photon energy is "+sum_particle.e()+" (GeV).");
            }
        }
		        	}
		        	reader.close();
		     	}
		      } 
	}
	void doDiHadrons(PhysicsEvent generic_Event)
	{
	  	for(int i=0;i<generic_Event.count();i++)
	  	{
	  		MyParticle part=(MyParticle)generic_Event.getParticle(i);
//            		System.out.println("matching mc particle index: " + part.matchingMCPartIndex);
                   		
   
	  		if(part.pid()==211)
	  		{ 			
	  				//	System.out.println("found proton");
	  			for(int j=0;j<generic_Event.count();j++)
	  			{
	  				MyParticle part2=(MyParticle)generic_Event.getParticle(j);
	  				//	Systefm.out.println("lookign at pid " + part2.pid());
	  				if(part2.pid()==(-211))
	  				{
	  					LorentzVector pionPair=new LorentzVector(part.px()+part2.px(),part.py()+part2.py(),part.pz()+part2.pz(),part.e()+part2.e());
	  					hDiPionMass.fill(pionPair.mass());
	  				    LorentzVector boostedPair=new LorentzVector(pionPair);	
	  				    boostedPair.boost(novel_fitter.gNBoost);
	  				    //we are in the breit frame now, so I guess that the transverse component
	  				    //should already be Pht
	  				    double px=boostedPair.vect().x();
	  				    double py=boostedPair.vect().y();
	  				    double phT= Math.sqrt(px*px+py*py);
	  				    hDiPionPPerp.fill(phT);
	  					
	  				}
	  			}
	  		}
	  	}
		
		
	}
	
	
	 public void plot() {
         EmbeddedCanvas can_lambda = new EmbeddedCanvas();
         can_lambda.setSize(1200,600);
         can_lambda.divide(2,2);
         can_lambda.setAxisTitleSize(24);
         can_lambda.setAxisFontSize(24);
         can_lambda.setTitleSize(24);
         can_lambda.cd(0);can_lambda.draw(hLambdaMass);
         can_lambda.cd(1);can_lambda.draw(hLambdaXf);
         can_lambda.cd(2);can_lambda.draw(hLambdaMassXfCut);
         can_lambda.cd(3);can_lambda.draw(hXfVsZ);
        // can_lambda.cd(1);can_e_ecal.draw(H_ESampl_ECal);
         //can_lambda.draw(g_m_ESampl_ECal,"same");
         //can_lambda.draw(g_s_ESampl_ECal,"same");
         can_lambda.save("lambda.png");
         
         EmbeddedCanvas can_piPi = new EmbeddedCanvas();
         can_piPi.setSize(1200,600);
         can_piPi.divide(2,2);
         can_piPi.setAxisTitleSize(24);
         can_piPi.setAxisFontSize(24);
         can_piPi.setTitleSize(24);
         can_piPi.cd(0);can_piPi.draw(hPiPiMass);
         can_piPi.cd(1);can_piPi.draw(hPiPiXf);
         can_piPi.cd(2);can_piPi.draw(hPiPiMassXfCut);
         can_piPi.cd(3);can_piPi.draw(hPiPiXfVsZ);
        // can_piPi.cd(1);can_e_ecal.draw(H_ESampl_ECal);
         //can_piPi.draw(g_m_ESampl_ECal,"same");
         //can_piPi.draw(g_s_ESampl_ECal,"same");
         can_piPi.save("piPiCandidates.png");
         
         
         EmbeddedCanvas can_Res = new EmbeddedCanvas();
         can_Res.setSize(1200,600);
         can_Res.divide(2,2);
         can_Res.setAxisTitleSize(24);
         can_Res.setAxisFontSize(24);
         can_Res.setTitleSize(24);
         can_Res.cd(0);can_Res.draw(hLambdaMassRes);
         can_Res.cd(1);can_Res.draw(hLambdaXfRes);
         can_Res.cd(2);can_Res.draw(hTrueLambdaMass);
         can_Res.cd(3);can_Res.draw(hTrueXf);
        // can_piPi.cd(1);can_e_ecal.draw(H_ESampl_ECal);
         //can_piPi.draw(g_m_ESampl_ECal,"same");
         //can_piPi.draw(g_s_ESampl_ECal,"same");
         can_Res.save("ResolutionStudies.png");
         
         EmbeddedCanvas can_dihad = new EmbeddedCanvas();
         can_dihad.setSize(1200,600);
         can_dihad.divide(2,2);
         can_dihad.setAxisTitleSize(24);
         can_dihad.setAxisFontSize(24);
         can_dihad.setTitleSize(24);
         can_dihad.cd(0);can_dihad.draw(hDiPionMass);
         can_dihad.cd(1);can_dihad.draw(hDiPionPPerp);
         can_dihad.save("dihadrons.png");
         
         
 }
	 protected void associateMCWithData(PhysicsEvent generic_Event, PhysicsEvent generic_EventMC)
	 {
		
		 for(int i=0;i<generic_Event.count();i++)
		 {	 
			 double minMomDiff=1000.0;
			 int minMomDiffIndex=-1;
			 MyParticle part=(MyParticle)generic_Event.getParticle(i);
			 double px=part.px();
			 double py=part.py();
			 double pz=part.pz();
			// System.out.println("data part px " + px + " py " + py + " pz " + pz);
			 
			 double mom = Math.sqrt(px*px+py*py+pz*pz);
             double theta = Math.toDegrees(Math.atan2(Math.sqrt(px*px+py*py), pz));
             double phi = Math.toDegrees(Math.atan2(py, px));
           //  System.out.println("data theta "+ theta + " phi " + phi + " mom " +mom);
			 for(int j=0;j<generic_EventMC.count();j++)
			 {	 				 
				 MyParticle partMC=(MyParticle)generic_EventMC.getParticle(j);
				 double pxMC=partMC.px();
				 double pyMC=partMC.py();
				 double pzMC=partMC.pz();
			
				 double momMC = Math.sqrt(pxMC*pxMC+pyMC*pyMC+pzMC*pzMC);
                 double thetaMC = Math.toDegrees(Math.atan2(Math.sqrt(pxMC*pxMC+pyMC*pyMC), pzMC));
                 double phiMC = Math.toDegrees(Math.atan2(pyMC, pxMC));
            	 if(i==0)
				 {
			//		 System.out.println("MC part px " + pxMC + " py " + pyMC + " pz " + pzMC);
			//		 System.out.println("mc teta "+ thetaMC + " phi " + phiMC + " mom " +momMC);
				 }
				 
      //    System.out.println("Looking at MC part "+ j + " rel momDiff: " + ((mom-momMC)/momMC));
                 if( Math.abs(mom-momMC) <0.025*momMC && Math.abs(theta-thetaMC)<1.0 && Math.abs(phi-phiMC)<5 )
                	 {
                		 if(Math.abs(mom-momMC)<minMomDiff)
                		 {
                			 minMomDiff=Math.abs(mom-momMC);	
                			 minMomDiffIndex=j;
                		 }
                	 }
			 }
			// System.out.println("associate mc part "+minMomDiffIndex + " with "+ i);
			 part.matchingMCPartIndex=minMomDiffIndex;
		 }
		
	
	 }

		 
	 
}
