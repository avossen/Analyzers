import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.io.*;
import java.nio.file.*;
import java.io.Serializable;
import java.util.*;
import java.util.StringTokenizer;

import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.groot.base.PadMargins;
import org.jlab.groot.base.Attributes;
import org.jlab.groot.base.GStyle;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
//import org.jlab.groot.math.Func2D;

import org.jlab.groot.data.GraphErrors;

import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.*;

public class PairReader {
static int debugEvent=51672129;
	protected int numGoodEvents;
	protected H1F phiRResolution;
	protected H1F thetaResolution;
	protected H1F zResolution;
	protected H1F hY;
	protected H1F hWy;
	protected H1F MResolution;

	protected H1F hPhiR;
	protected H1F hPhiH;
	protected H1F hG1PAngle;

	protected H2F hPhiRVsPhiH;
	protected H2F hPhiRVsG1PAngle;

	// protected H1F hPhiR;

	protected H1F hZ;
	protected H1F hXf;

	protected H1F hM;

	protected H1F hAsymDistVsRun;
	protected H1F hAsymDistVsRunG1P;
	protected H2F phiRVsZResolution;
	protected H2F MVsZResolution;
	protected H2F hQ2VsX;

	// is the halfway plate in
	protected boolean hwpIn = false;
	static boolean isMC = false;
	
	// static boolean isMC = true;
	//static boolean isHrayr = true;
	static boolean isHrayr = false;
	static int numPhiBins = 8;
	static int numPhiBins2D = 8;
	// this is driven by the run numbers
	static int numOutbendingRuns=68;
	static int maxKinBins = 42+numOutbendingRuns;
	int numPairs=0;
	// arrays for asymmetry computation. Let's just to pi+pi for now
	// so this is indexed in the kinBin, spin state, phi bin
	protected float[][][][] counts;
	protected float[][][][] countsG1P;

	protected float[][][][][] counts_2D;
	protected float[][][][][] countsG1P_2D;

	protected double[][] m_pairsWithMatch;
	protected double[][] m_pairsWithP1Match;
	protected double[][] m_pairsWithP2Match;
	protected double[][] m_pairsWithSingleMatch;
	protected double[][] m_pairsWOMatch;

	protected double[][] meanKin;
	// to calculate multiplicity per event
	protected double[] numEvts;
	// to compute kin variables vs run number (and maybe others)
	protected double[][] meanQ2;
	protected double[][] meanPt;
	protected double[][] meanM;
	protected double[][] meanHadMom1;
	protected double[][] meanHadMom2;
	protected double[][] meanHadPhi1;
	protected double[][] meanHadPhi2;
	protected double[][] meanHadTheta1;
	protected double[][] meanHadTheta2;
	protected double[] meanHadPairMult;

	protected H1F theta1;
	protected H1F theta2;
	protected H1F phi1;
	protected H1F phi2;

	protected float[][] meanWy;
	protected float[][][] kinCount;

	protected ArrayList<Double> phiBins;
	protected ArrayList<Double> phiBins2D;

	public static void main(String[] args) {
		GStyle.getGraphErrorsAttributes().setMarkerStyle(0);
		GStyle.getGraphErrorsAttributes().setMarkerColor(3);
		GStyle.getGraphErrorsAttributes().setMarkerSize(7);
		GStyle.getGraphErrorsAttributes().setLineColor(3);
		GStyle.getGraphErrorsAttributes().setLineWidth(3);
		GStyle.getAxisAttributesX().setTitleFontSize(34);
		GStyle.getAxisAttributesX().setLabelFontSize(25);
		GStyle.getAxisAttributesY().setTitleFontSize(32);
		GStyle.getAxisAttributesY().setLabelFontSize(25);

		PairReader myReader = new PairReader();
		myReader.initialize();
		myReader.analyze(args);
		myReader.savePlots();

	}

	public void initialize() {
		theta1 = new H1F("theta1", "theta1", 100, 0, 1.5);
		theta2 = new H1F("theta2", "theta2", 100, 0, 1.5);

		phi1 = new H1F("phi1", "phi1", 100, -Math.PI, Math.PI);
		phi2 = new H1F("phi2", "phi2", 100, -Math.PI, Math.PI);

		phiRResolution = new H1F("phiRResolution", "phiRResolution", 100, -0.3, 0.3);
		thetaResolution = new H1F("thetaResolution", "thetaResolution", 100, -0.3, 0.3);
		zResolution = new H1F("zResolution", "zResolution", 100, -0.3, 0.3);
		MResolution = new H1F("MResolution", "MResolution", 100, -0.3, 0.3);
		hY = new H1F("Y", "Y", 100, -0.3, 1.2);
		hWy = new H1F("Wy", "Wy", 100, -0.3, 1.2);
		hZ = new H1F("Z", "Z", 100, 0, 1.0);
		hM = new H1F("M", "M", 100, 0, 3.0);
		hAsymDistVsRun = new H1F("asymDistVsRun", "asymDistVsRun", 20, -5, 5);
		hAsymDistVsRunG1P = new H1F("asymDistVsRunG1P", "asymDistVsRunG1P", 20, -5, 5);
		hXf = new H1F("xF", "xF", 100, -1.0, 1.0);
		hPhiR = new H1F("phiR", "phiR", 100, 0, 2 * Math.PI);
		hPhiH = new H1F("phiH", "phiH", 100, 0, 2 * Math.PI);
		hG1PAngle = new H1F("g1pAngle", "g1pAngle", 100, 0, 2 * Math.PI);

		this.hPhiRVsG1PAngle = new H2F("phiRvsG1PAngle", "phiRVsG1PAngle", 30, 0, 2 * Math.PI, 30, 0, 2 * Math.PI);
		this.hPhiRVsPhiH = new H2F("phiRvsPhiH", "phiRVsPhiH", 30, 0, 2 * Math.PI, 30, 0, 2 * Math.PI);

		MVsZResolution = new H2F("MVsZResolution", "MVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		phiRVsZResolution = new H2F("phiRVsZResolution", "phiRVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		hQ2VsX = new H2F("Q2VsX", "Q2VsX", 20, 0.0, 1.0, 20, 0.0, 12);
		phiBins = new ArrayList<Double>();
		phiBins2D = new ArrayList<Double>();
		counts = new float[Binning.numKinBins][2][maxKinBins][numPhiBins];
		countsG1P = new float[Binning.numKinBins][2][maxKinBins][numPhiBins];
		
		counts_2D = new float[Binning.numKinBins][2][maxKinBins][numPhiBins2D][numPhiBins2D];
		countsG1P_2D = new float[Binning.numKinBins][2][maxKinBins][numPhiBins2D][numPhiBins2D];

		meanKin = new double[Binning.numKinBins][maxKinBins];
		m_pairsWithMatch = new double[Binning.numKinBins][maxKinBins];
		m_pairsWithP1Match = new double[Binning.numKinBins][maxKinBins];
		m_pairsWithP2Match = new double[Binning.numKinBins][maxKinBins];
		m_pairsWithSingleMatch = new double[Binning.numKinBins][maxKinBins];
		m_pairsWOMatch = new double[Binning.numKinBins][maxKinBins];
		meanQ2 = new double[Binning.numKinBins][maxKinBins];
		meanHadPairMult = new double[maxKinBins];
		numEvts = new double[maxKinBins];
		meanHadMom1 = new double[Binning.numKinBins][maxKinBins];
		meanHadMom2 = new double[Binning.numKinBins][maxKinBins];
		meanHadTheta1 = new double[Binning.numKinBins][maxKinBins];
		meanHadTheta2 = new double[Binning.numKinBins][maxKinBins];
		meanHadPhi1 = new double[Binning.numKinBins][maxKinBins];
		meanHadPhi2 = new double[Binning.numKinBins][maxKinBins];
		meanPt = new double[Binning.numKinBins][maxKinBins];
		meanM = new double[Binning.numKinBins][maxKinBins];

		meanWy = new float[Binning.numKinBins][maxKinBins];
		// also needed for relative luminosity
		kinCount = new float[Binning.numKinBins][2][maxKinBins];

		// can't do this with multidimensional arrays (but default is already 0)
		// Arrays.fill(counts, (float) 0.0);
		// Arrays.fill(meanKin, (float) 0.0);
		// Arrays.fill(kinCount, (float) 0.0);

		for (int i = 0; i < numPhiBins; i++) {
			System.out.println("adding phi bin " + ((i + 1) * 2 * Math.PI / numPhiBins));
			phiBins.add(((i + 1) * 2 * Math.PI / numPhiBins));
		}
		for (int i = 0; i < numPhiBins2D; i++) {
			System.out.println("adding phi bin " + ((i + 1) * 2 * Math.PI / numPhiBins2D));
			phiBins2D.add(((i + 1) * 2 * Math.PI / numPhiBins2D));
		}

	}

	public void savePlots() {
		EmbeddedCanvas can_piPi = new EmbeddedCanvas();
		can_piPi.setSize(800, 600);
		can_piPi.divide(2, 2);
		// can_piPi.setAxisTitleSize(24);
		// can_piPi.setAxisFontSize(24);

		// can_piPi.setTitleSize(24);
		can_piPi.cd(0);

		phiRResolution.setTitleX("phiR Resolution");
		zResolution.setTitleX("z Resolution");
		thetaResolution.setTitleX("theta Resolution");
		MResolution.setTitleX("M Resolution");
		// can_piPi.getPad(0).getAxisX().setTitle("phiR Resolution");
		// can_piPi.getPad(1).getAxisX().setTitle("theta Resolution");
		// can_piPi.getPad(2).getAxisX().setTitle("z Resolution");
		// can_piPi.getPad(3).getAxisX().setTitle("M Resolution");
		can_piPi.draw(phiRResolution);
		can_piPi.cd(1);
		;
		can_piPi.draw(thetaResolution);
		can_piPi.cd(2);
		;
		can_piPi.draw(zResolution);
		can_piPi.cd(3);
		;
		can_piPi.draw(MResolution);
		can_piPi.save("resolutions.png");

		EmbeddedCanvas canTP = new EmbeddedCanvas();
		canTP.setSize(800, 600);
		canTP.divide(2, 2);
		canTP.cd(0);
		canTP.getPad(0).getAxisY().setTitle("Theta1");
		theta1.divide((double) numGoodEvents);
		canTP.draw(theta1);

		canTP.cd(1);
		canTP.getPad(1).getAxisY().setTitle("Theta1");
		theta2.divide((double) numGoodEvents);
		canTP.draw(theta2);

		canTP.cd(2);
		canTP.getPad(2).getAxisY().setTitle("phi1");
		phi1.divide((double) numGoodEvents);
		canTP.draw(phi1);

		canTP.cd(3);
		canTP.getPad(3).getAxisY().setTitle("phi2");
		phi2.divide((double) numGoodEvents);
		canTP.draw(phi2);

		canTP.save("thetaPhiDist.png");

		EmbeddedCanvas can2D = new EmbeddedCanvas();
		can2D.setSize(800, 600);
		can2D.divide(2, 1);
		can2D.cd(0);
		can_piPi.getPad(0).getAxisY().setTitle("phiR Resolution");
		can_piPi.getPad(0).getAxisX().setTitle("z");
		can_piPi.getPad(1).getAxisY().setTitle("M Resolution");
		can_piPi.getPad(0).getAxisX().setTitle("z");
		can2D.draw(this.phiRVsZResolution);
		can2D.cd(1);
		can2D.draw(this.MVsZResolution);
		can2D.save("twoDResolutions.png");

		EmbeddedCanvas canAngles = new EmbeddedCanvas();
		canAngles.setSize(800, 600);
		canAngles.divide(3, 2);
		canAngles.cd(0);
		canAngles.getPad(0).getAxisX().setTitle("phiR");
		canAngles.draw(this.hPhiR);

		canAngles.cd(1);
		canAngles.getPad(1).getAxisX().setTitle("phiH");
		canAngles.draw(this.hPhiH);

		canAngles.cd(2);
		canAngles.getPad(2).getAxisX().setTitle("phiH-PhiH");
		canAngles.draw(this.hG1PAngle);

		canAngles.cd(3);
		canAngles.getPad(3).getAxisX().setTitle("phiRVsG1PAngle");
		canAngles.draw(this.hPhiRVsG1PAngle);

		canAngles.cd(4);
		canAngles.getPad(4).getAxisX().setTitle("phiRVsPhiH");
		canAngles.draw(this.hPhiRVsPhiH);

		// canAngles.getPad(1).getAxisX().setTitle("phiH");
		// canAngles.cd(1);
		// canAngles.draw(this.hPhiH);
		canAngles.save("angles.png");

		EmbeddedCanvas canKin = new EmbeddedCanvas();
		canKin.setSize(800, 600);
		canKin.divide(2, 3);
		canKin.cd(0);
		canKin.getPad(0).getAxisX().setTitle("z");
		canKin.draw(this.hZ);
		canKin.cd(1);
		canKin.getPad(1).getAxisX().setTitle("M");

		canKin.draw(this.hM);
		canKin.getPad(2).getAxisX().setTitle("xF");
		canKin.cd(2);
		canKin.draw(this.hXf);
		canKin.cd(3);
		canKin.getPad(3).getAxisX().setTitle("x");
		canKin.getPad(3).getAxisY().setTitle("Q2");
		canKin.getPad(3).getAxisZ().setLog(true);
		canKin.draw(this.hQ2VsX);
		canKin.cd(4);
		canKin.getPad(4).getAxisX().setTitle("Y");
		canKin.draw(this.hY);
		canKin.cd(5);
		canKin.getPad(5).getAxisX().setTitle("Wy");
		canKin.draw(this.hWy);
		canKin.save("diHadKins.png");

		// fit with gaussian
		F1D fG = new F1D("fG", "[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
		// Func2D f2G = new Func2D( ("fG","[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
		fG.setParameter(0, 5.0);
		fG.setParameter(1, 0.0);
		fG.setParameter(2, 1.0);

		F1D fG2 = new F1D("fG2", "[amp]*gaus(x,[mean],[sigma])", -5.0, 5.0);
		fG2.setParameter(0, 5.0);
		fG2.setParameter(1, 0.0);
		fG2.setParameter(2, 1.0);

		// DataFitter2D.fit(Func2D,H2F);
		DataFitter.fit(fG, this.hAsymDistVsRun, "Q");
		EmbeddedCanvas canStability = new EmbeddedCanvas();
		canStability.setSize(800, 600);
		canStability.divide(2, 1);
		canStability.cd(0);
		canStability.getPad(0).getAxisX().setTitle("normalized asymmetry");
		canStability.draw(this.hAsymDistVsRun);
		canStability.draw(fG, "same");

		System.out.println("Gaus fit to asyms mean: " + fG.parameter(1).value() + " +- " + fG.parameter(1).error()
				+ " sigma: " + fG.parameter(2).value() + " +- " + fG.parameter(2).error());
		System.out.println("chi2/ndf: " + fG.getChiSquare() / fG.getNDF());

		DataFitter.fit(fG2, this.hAsymDistVsRunG1P, "Q");
		canStability.cd(1);
		canStability.getPad(1).getAxisX().setTitle("normalized asymmetry G1P");
		canStability.draw(this.hAsymDistVsRunG1P);
		canStability.draw(fG2, "same");
		System.out.println("Gaus fit to asyms G1P mean: " + fG2.parameter(1).value() + " +- " + fG2.parameter(1).error()
				+ " sigma: " + fG2.parameter(2).value() + " +- " + fG2.parameter(2).error());
		System.out.println("chi2/ndf: " + fG2.getChiSquare() / fG2.getNDF());
		canStability.save("asymStability.png");

	}

	public void analyze(String[] args) {
		int pairsWithMatch = 0;
		int pairsWOMatch = 0;
		numGoodEvents = 0;
		// TODO Auto-generated method stub
		AsymData m_asymData = null;
		File folder = new File(args[0]);
		File[] listOfFiles = folder.listFiles();
		for (int iF = 0; iF < listOfFiles.length; iF++) {
			if (listOfFiles[iF].isFile()) {
				System.out.println("File " + listOfFiles[iF].getName());
				PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:*.{srn}");

				Path filename = Paths.get(listOfFiles[iF].getName());
				if (matcher.matches(filename)) {
					//System.out.println("looing at " + filename);
					// Deserialization
					// check run number from filename, so we can determine halfway plate
					// position
					// is also in the dst's but not in my srn files (yet)
					int runNumber = 0;
					boolean foundRunNumber=false;
					StringTokenizer st = new StringTokenizer(listOfFiles[iF].getName(), "_");
					for (int i = 0; st.hasMoreTokens(); i++) {
						//let's see if the run number is the first thing like in the
						//new skim files	
						String t1 = st.nextToken();
						try
						{
							runNumber = Integer.parseInt(t1);		
							if(runNumber>999 && runNumber<10000)
								break;
						}
						catch(NumberFormatException ex)
						{
							runNumber=0;
						}
						// System.out.println("token " + i + ": "+ t1);
						if (t1.contentEquals("out")) {
							// System.out.println("passed first check");
							String t2 = st.nextToken();
							// System.out.println("next token " + i + ": "+ t2);
							if (t2.contentEquals("clas")) {

								String t3 = st.nextToken(".");
								// need to strip off leading "_" if we use "." token.
								t3 = t3.substring(1, t3.length());
								// System.out.println("getting run number "+t3);
								runNumber = Integer.parseInt(t3);

							}
						}
					}
					
					try {
						// Reading the object from a file
						ObjectInputStream in = null;
						FileInputStream file = null;
						try {
							file = new FileInputStream(args[0] + "/" + filename.toString());
							in = new ObjectInputStream(file);
						} catch (IOException e) {
							System.out.println("Caught IO exception: " + e.getMessage());
							throw e;
						}
						// Method for deserialization of object
						m_asymData = (AsymData) in.readObject();
						System.out.println("got " + m_asymData.eventData.size() + " events");
						for (EventData evtData : m_asymData.eventData) {
							// helicity has to be fixed per event
							double rndSpin = 0.0;
							runNumber=evtData.runNr;
							//System.out.println("looking at run: " + runNumber);
							if (!isMC) {
								try {
									hwpIn = getHWPInfo(runNumber);
								} catch (Exception e) {

									System.out.println("wrong runnumber " + runNumber);
								}
							}
							
							
							//these seem to be outliers using 700 MeV as the cutoff for outliers
							if(runNumber==4308||runNumber==4302 || runNumber==4032|| runNumber ==4013 || runNumber ==4016|| runNumber==4017|| runNumber==4069||runNumber==4071||runNumber==4305||runNumber==4015	)
							{
								//continue;
							}
							//this run seems to have bad chi2 (2.2)
							if(runNumber==4071)
							{
							//	continue;
							}
							//use only runs marked as production runs
							if(!isMC)
							{
							if(runNumber>4115||runNumber<3910)
							{
								continue;
							}
							}
								//these seem to be bad with the 1 GeV cut
							if(runNumber==4069 || runNumber==4305 || runNumber==4071 || runNumber==4069|| runNumber==4039|| runNumber==4027||runNumber==4038|| runNumber==4302|| runNumber==4032||runNumber==4013)
							{
								//continue;
							}
							if (runNumber > 4200 && !isMC) {
								// continue;
							}
							if (runNumber < 3000 && !isMC) {
								// continue;
							}
							if ((runNumber == 4020 || runNumber == 4301) && !isMC) {
								// continue;
							}
							// so that we don't remove that run when plotting

							//System.out.println("run nr: " + runNumber);
							if (isMC) {
								// even MC has helicity-->only the stuff from Hrayr, the other one not,
								// so for now now need to set from random
								
									rndSpin = Math.random() * 2.0 - 1.0;
							}
						
							//System.out.println("beam hlicity: " + evtData.beamHelicity);
							int helicityIndex = 0;
							if (!isMC || isHrayr)
							{
//System.out.println("beam helicity: " + evtData.beamHelicity);
								if (evtData.beamHelicity > 0) {
									helicityIndex = 1;
									if(!isHrayr)
									{
										rndSpin = 1.0;
									}
								} else {
									if(!isHrayr)
									{
										rndSpin = -1.0;
									}
								}
							}
							if (isMC && runNumber == 0) {
								// let's have two runnumbers, otherwise groot doesn't make the graphs
								if (rndSpin > 0) {
									runNumber = 1;
								} else {
									// second bin
									runNumber = 4014;
								}

							}
						//	System.out.println("looking at run number: " + runNumber);
							if (isMC && !isHrayr) {
								if (rndSpin < 0) {
									helicityIndex = 1;
								} else {
									helicityIndex = 0;
								}
								
							}
							// flip helicities
							if (hwpIn) {
								System.out.println("change helicyt index");
								if (helicityIndex == 1) {
									helicityIndex = 0;
									rndSpin = -1.0;
								} else {
									helicityIndex = 1;
									rndSpin = 1.0;
								}
							}
							//System.out.println("helicity: " + helicityIndex);
							hQ2VsX.fill(evtData.x, evtData.Q2);
							// kinematic factor W(y) (asymmetry is ~W(y) x e(x)
							double Wy = 0;
							if (evtData.y > 0)
								Wy = 2 * evtData.y * Math.sqrt(1 - evtData.y);
							else
								continue;

							hY.fill(evtData.y);
							hWy.fill(Wy);
							double nu = evtData.y * 10.6;
							int runBin = Binning.RunBinning.getBin(0.0, 0.0, evtData.x, runNumber);
							meanHadPairMult[runBin] += evtData.pairData.size();
							numEvts[runBin] += 1;

							numGoodEvents++;
							if(evtData.eventNr==debugEvent)
							{
								System.out.println("found event!!\n");
							}
							
							
							boolean foundGoodPair=false;
							for (HadronPairData pairData : evtData.pairData) {
								theta1.fill(pairData.theta1);
								theta2.fill(pairData.theta2);

								phi1.fill(pairData.phi1);
								phi2.fill(pairData.phi2);
								if (pairData.theta1 < 0.17 || pairData.theta2 < 0.17) {
								//if (pairData.theta1 < 0.087 || pairData.theta2 < 0.087) {
									continue;
								}
								// acceptance of forward detector
								if (pairData.theta1 > 0.5236 || pairData.theta2 > 0.5236) {
									//continue;
								}
								
								//float E1=Math.sqrt((double)(pairData.mom1*pairData.mom1+0.14*0.14));
								//float E2=Math.sqrt((double)(pairData.mom2*pairData.mom2+0.14*0.14));
								//System.out.println("found track with mom: "+pairData.mom1 + " "+ pairData.mom2);
								if (pairData.mom1 <1.0|| pairData.mom2 < 1.0) {
									//System.out.println("found track with mom: "+pairData.mom1 + " "+ pairData.mom2);
									continue;
								}
								if(pairData.M< 0.8)	
								{
								 //continue;
								}
								hXf.fill(pairData.xF);
								hZ.fill(pairData.z);
								hM.fill(pairData.M);
								// System.out.println("helicity: " + evtData.beamHelicity);
								if (pairData.z < 0.2 || pairData.xF < 0)
								{
									continue;
								}

								double eh1 = Math.sqrt(pairData.mom1 * pairData.mom1 + 0.14 * 0.14);
								double eh2 = Math.sqrt(pairData.mom2 * pairData.mom2 + 0.14 * 0.14);
								double z1 = eh1 / nu;
								double z2 = eh2 / nu;
								if (z1 < 0.1 || z2 < 0.1)
								{
									continue;
								}
								///Debug print
								System.out.print("z1: " + z1 + " z2 "+z2);
								System.out.println(" pion1 mom: "+ pairData.mom1+" pion2 mom: " +pairData.mom2);
									
						
								///--->end
								// if(pairData.M<0.7)
								// continue;
								numPairs++;
								double weight = 1.0;
								foundGoodPair=true;
								if (pairData.hasMC) {
									weight = getWeight(pairData.matchingMCPair, rndSpin);
								}
								int phiBin = Binning.getBin(phiBins, pairData.phiR);
								int phiBin2D = Binning.getBin(phiBins2D, pairData.phiR);
								double g1Pangle = pairData.phiH - pairData.phiR;
								if (g1Pangle < 0)
									g1Pangle = g1Pangle + 2 * Math.PI;
								hPhiR.fill(pairData.phiR);
								hPhiH.fill(pairData.phiH);
								hG1PAngle.fill(g1Pangle);

								this.hPhiRVsG1PAngle.fill(g1Pangle, pairData.phiR);
								this.hPhiRVsPhiH.fill(pairData.phiH, pairData.phiR);
								int phiBinG1P = Binning.getBin(phiBins, g1Pangle);
								int phiHBin2D = Binning.getBin(phiBins2D, pairData.phiH);
								// System.out.println("phiR: "+pairData.phiR+ " bin: "+ phiBin);
								for (Binning binningType : EnumSet.allOf(Binning.class)) {
									int iBin = binningType.getBin(pairData.M, pairData.z, evtData.x, runNumber);
									// int phiBin = binningType.getBin(phiBins, pairData.phiR);

									// System.out.println("kin bin " + binningType.name() + " m: " +pairData.M +" z:
									// "+ pairData.z +" x: " + evtData.x+ " bin: "+ iBin);

									if (iBin < 0) {
										System.out.println("kinematic bin too small " + binningType.name() + " m: "
												+ pairData.M + " z: " + pairData.z + " x: " + evtData.x);
									}
									if (phiBin < 0) {
										System.out.println("phi value " + pairData.phiR + " out of bounds");
									}

									// do it that way since e.g. for mc data we get -99

									counts[binningType.binType][helicityIndex][iBin][phiBin] += weight;
									countsG1P[binningType.binType][helicityIndex][iBin][phiBinG1P] += weight
											* pairData.pTBreit / pairData.M;

									counts_2D[binningType.binType][helicityIndex][iBin][phiBin2D][phiHBin2D] += weight;
									countsG1P_2D[binningType.binType][helicityIndex][iBin][phiBin2D][phiHBin2D] += weight
											* pairData.pTBreit / pairData.M;

									kinCount[binningType.binType][helicityIndex][iBin] += weight;
									meanWy[binningType.binType][iBin] += Wy * weight;
									meanQ2[binningType.binType][iBin] += evtData.Q2 * weight;
									meanPt[binningType.binType][iBin] += pairData.pTBreit * weight;
									meanM[binningType.binType][iBin] += pairData.M * weight;
									meanHadPhi1[binningType.binType][iBin] += pairData.phi1 * weight;
									meanHadPhi2[binningType.binType][iBin] += pairData.phi2 * weight;

									meanHadTheta1[binningType.binType][iBin] += pairData.theta1 * weight;
									meanHadTheta2[binningType.binType][iBin] += pairData.theta2 * weight;

									meanHadMom1[binningType.binType][iBin] += pairData.mom1 * weight;
									meanHadMom2[binningType.binType][iBin] += pairData.mom2 * weight;

									if (pairData.hasMC) {
										m_pairsWithMatch[binningType.binType][iBin] += weight;
									} else {
										m_pairsWOMatch[binningType.binType][iBin] += weight;
									}

									if (binningType == Binning.MBinning) {
										meanKin[binningType.binType][iBin] += pairData.M * weight;

									}
									if (binningType == Binning.ZBinning) {
										meanKin[binningType.binType][iBin] += pairData.z * weight;

									}
									if (binningType == Binning.XBinning)
										meanKin[binningType.binType][iBin] += evtData.x * weight;
									if (binningType == Binning.RunBinning) {
										meanKin[binningType.binType][iBin] += runNumber * weight;

									}
								}

								// hPhiH.fill(pairData.ph);
								if (pairData.hasMC) {

									phiRResolution.fill((pairData.phiR - pairData.matchingMCPair.phiR) / pairData.phiR);
									pairsWithMatch++;

									// System.out.println("pair theta:" + pairData.theta +" matching theta : " +
									// pairData.matchingMCPair.theta );
									thetaResolution
											.fill((pairData.theta - pairData.matchingMCPair.theta) / pairData.theta);
									zResolution.fill((pairData.z - pairData.matchingMCPair.z) / pairData.z);
									MResolution.fill((pairData.M - pairData.matchingMCPair.M) / pairData.M);
									MVsZResolution.fill(pairData.z,
											(pairData.M - pairData.matchingMCPair.M) / pairData.M);
									phiRVsZResolution.fill(pairData.z,
											(pairData.phiR - pairData.matchingMCPair.phiR) / pairData.phiR);

									System.out.println("Hadron pair with phiR: " + pairData.phiR + " has mc partner: "
											+ pairData.matchingMCPair.phiR);
								} else {

									pairsWOMatch++;
								}
								// System.out.println("looking at hadron pair with z: "+ pairData.z);
							}
							if(foundGoodPair)
							{
								System.out.print("event nr " + evtData.eventNr+ " runnr "+evtData.runNr);
								System.out.println(" Q2 " + evtData.Q2+ " y "+evtData.y +" x "+evtData.x);
								System.out.println("\n");
							}
						}
						
						in.close();
						file.close();

						System.out.println("Object has been deserialized ");
					} catch (IOException ex) {
						System.out.println("IOException is caught");
					} catch (ClassNotFoundException ex) {
						System.out.println("ClassNotFoundException is caught");
					}

				}
			}
		}
		// ran over all files, now fit
		dump2DToFile(false);
		dump2DToFile(true);
		doFits(false);
		doFits(true);
		System.out.println("Pairs with match " + pairsWithMatch + " pairs without: " + pairsWOMatch + " percentage: "
				+ pairsWithMatch / (float) (pairsWithMatch + pairsWOMatch));
		System.out.println("num pairs: "+numPairs);
	}

	
	void writeH2ToFile(H2F h2)
	{
		String filename=h2.getName()+".txt";
	
		Writer writer=null;	
		try 
		{
		writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "utf-8")); 
		writer.write(h2.getXAxis().getNBins()+" "+h2.getYAxis().getNBins()+" ");
		writer.write(h2.getXAxis().min()+" "+h2.getYAxis().min()+" ");
		writer.write(h2.getXAxis().max()+" "+h2.getYAxis().max()+" ");
		for(int i=0;i<h2.getXAxis().getNBins();i++)
		{
			for(int j=0;j<h2.getXAxis().getNBins();j++)
			{
				double val=h2.getBinContent(i,j);
				writer.write(Float.toString((float)val)+" ");	
			}
		}
		writer.write("\n");
		}
		catch(IOException ex)
		{
			
		}
		finally
		{
			try {writer.close();} catch (Exception ex) {/*ignore*/}
		}
		
		
	}
	void writeH1ToFile(H1F h1)
	{
		String filename=h1.getName()+".txt";
		Writer writer=null;
		try 
		{
		writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "utf-8")); 
		writer.write(h1.getXaxis().getNBins()+" ");
		writer.write(h1.getXaxis().min()+" "+h1.getXaxis().max()+" ");
		for(int i=0;i<h1.getXaxis().getNBins();i++)
		{
			double val=h1.getBinContent(i);
			writer.write(Float.toString((float)val)+" ");	
		}
		writer.write("\n");
		}
		catch(IOException ex)
		{
			
		}
		finally
		{
			try {writer.close();} catch (Exception ex) {/*ignore*/}
		}
	}
	
	void dump2DToFile(boolean isG1P)
	{
		String filename;
		Writer writer=null;
		float loc_counts[][][][][];
//		// g1p counts are weighted
		if (isG1P)
			loc_counts = countsG1P_2D;
		else
			loc_counts = counts_2D;

		if(isG1P)
		{
			filename="countsOutG1P.txt";
		}
		else
		{
			filename="countsOut.txt";
		}
		try 
		{
			writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "utf-8")); 			
			for (Binning binningType : EnumSet.allOf(Binning.class)) {
				writer.write(binningType.getBinningName()+" "+Integer.toString(binningType.getNumBins())+" "+Integer.toString(numPhiBins)+" ");
				//dump mean values for each kin bin
				for (int iKinBin = 0; iKinBin < binningType.getNumBins(); iKinBin++) {
					double normCount=(kinCount[binningType.getBinType()][0][iKinBin]
						+ kinCount[binningType.getBinType()][1][iKinBin]);
					double xVal = this.meanKin[binningType.binType][iKinBin]/normCount;	
					writer.write(Float.toString((float)xVal)+" ");
					double wyFactor = this.meanWy[binningType.binType][iKinBin]
							/ normCount;
					writer.write(Float.toString((float)wyFactor)+" ");
				}
				writer.write("\n");
				for (int iKinBin = 0; iKinBin < binningType.getNumBins(); iKinBin++) {
					float r=0;
					if (kinCount[binningType.getBinType()][1][iKinBin] > 0) {
						r = kinCount[binningType.getBinType()][0][iKinBin]
								/ kinCount[binningType.getBinType()][1][iKinBin];
					} else {
						System.out.println("dumpfile no counts (r) for kin bin " + iKinBin);
					}
					writer.write(Float.toString(r)+" ");
					for (int iAngBin = 0; iAngBin < numPhiBins; iAngBin++) {
						for (int iAngBin2 = 0; iAngBin2 < numPhiBins; iAngBin2++) {	
							for(int ipol=0;ipol<2;ipol++)
							{
								writer.write(Float.toString(loc_counts[binningType.getBinType()][ipol][iKinBin][iAngBin][iAngBin2]));
								writer.write(" ");
							}
	
						}
					}
				}
				writer.write("\n");
			}

		}
			catch(IOException ex)
			{
				
			}
			finally
			{
				try {writer.close();} catch (Exception ex) {/*ignore*/}
			}
		
	}
	
	void doFits2D(boolean doG1P) {
//		My2DFit f2 = new My2DFit("myFit", 2, 0, 2 * Math.PI, 0, 2 * Math.PI);
//
//		float loc_counts[][][][][];
//		// g1p counts are weighted
//		if (doG1P)
//			loc_counts = countsG1P_2D;
//		else
//			loc_counts = counts_2D;
//
//	
//
//		for (Binning binningType : EnumSet.allOf(Binning.class)) {
//
//			double vals[] = new double[binningType.getNumBins()];
//			double m_normCount[] = new double[binningType.getNumBins()];
//			double valErrs[] = new double[binningType.getNumBins()];
//			double xVals[] = new double[binningType.getNumBins()];
//			double xErrVals[] = new double[binningType.getNumBins()];
//			String baseTitle;
//			if (doG1P)
//				baseTitle = "G1P_";
//			else
//				baseTitle = "ex_";
//			for (int iKinBin = 0; iKinBin < binningType.getNumBins(); iKinBin++) {
//				f2.parameter(0).setValue(0.0);
//				f2.parameter(1).setValue(0.0);
//				
//				
//				String s = "myFitGraph_" + baseTitle + binningType.getBinningName() + "_bin" + iKinBin;
//				//the histogram is twice as wide in x direction to do the simultaneous fit
//				H2F hFit=new H2F("s","s",2*numPhiBins,0.0,4*Math.PI,numPhiBins,0,2*Math.PI);
//				
//
//				for (int iAngBin = 0; iAngBin < numPhiBins; iAngBin++) {
//					for (int iAngBin2 = 0; iAngBin2 < numPhiBins; iAngBin2++) {
//						double x = 0;
//						double y = 0;
//						double ex = 0;
//						double ey = 0;
//
//						double r = 1.0;
//						if (kinCount[binningType.getBinType()][1][iKinBin] > 0) {
//							r = kinCount[binningType.getBinType()][0][iKinBin]
//									/ kinCount[binningType.getBinType()][1][iKinBin];
//						} else {
//							System.out.println("no counts (r) for phi bin " + iAngBin + " " + s);
//						}
//						double N1 = loc_counts[binningType.getBinType()][0][iKinBin][iAngBin][iAngBin2];
//						double N2 = loc_counts[binningType.getBinType()][1][iKinBin][iAngBin][iAngBin2];
//						
//						
//					}
//					DataFitter.fit(f1, g, "Q");
//					vals[iKinBin] = f1.parameter(0).value();
//					valErrs[iKinBin] = f1.parameter(0).error();
//
//					if ((kinCount[binningType.getBinType()][0][iKinBin]
//							+ kinCount[binningType.getBinType()][1][iKinBin]) > 0) {
//						double normCount = (kinCount[binningType.getBinType()][0][iKinBin]
//								+ kinCount[binningType.getBinType()][1][iKinBin]);
//						xVals[iKinBin] = this.meanKin[binningType.binType][iKinBin] / normCount;
//						m_normCount[iKinBin] = normCount;
//						if (binningType == Binning.RunBinning) {
//							System.out.print("run bin" + iKinBin + ", meanKin: ");
//							System.out.print(this.meanKin[binningType.binType][iKinBin]);
//							System.out.println(", normFactor: " + normCount + " xval: " + xVals[iKinBin]);
//						}
//
//						double wyFactor = this.meanWy[binningType.binType][iKinBin] / normCount;
//					}
//				}
//			}
//		}
	}

	void doFits(boolean doG1P) {		
		float loc_counts[][][][];
		if (doG1P)
			loc_counts = countsG1P;
		else
			loc_counts = counts;
		// still no depolarisation factor and beam polarization
		// same function for both, e(x) and G1P
		F1D f1 = new F1D("f1", "0+[amp]*sin(x)", 0.0, 2 * Math.PI);
		f1.setParameter(0, 0.0);
		// f1.setParameter(1, 0.0);
		for (Binning binningType : EnumSet.allOf(Binning.class)) {

			double vals[] = new double[binningType.getNumBins()];
			double m_normCount[] = new double[binningType.getNumBins()];
			double valErrs[] = new double[binningType.getNumBins()];
			double xVals[] = new double[binningType.getNumBins()];
			double xErrVals[] = new double[binningType.getNumBins()];
			String baseTitle;
			if (doG1P)
				baseTitle = "G1P_";
			else
				baseTitle = "ex_";
			for (int iKinBin = 0; iKinBin < binningType.getNumBins(); iKinBin++) {
				String s = "myFitGraph_" + baseTitle + binningType.getBinningName() + "_bin" + iKinBin;
				GraphErrors g = new GraphErrors(s);

				for (int iAngBin = 0; iAngBin < numPhiBins; iAngBin++) {
					double x = 0;
					double y = 0;
					double ex = 0;
					double ey = 0;

					double r = 1.0;
					if (kinCount[binningType.getBinType()][1][iKinBin] > 0) {
						r = kinCount[binningType.getBinType()][0][iKinBin]
								/ kinCount[binningType.getBinType()][1][iKinBin];
					} else {
						System.out.println("no counts (r) for phi bin " + iAngBin + " " + s);
					}
					double N1 = loc_counts[binningType.getBinType()][0][iKinBin][iAngBin];
					double N2 = loc_counts[binningType.getBinType()][1][iKinBin][iAngBin];
					if ((N1 + r * N2) > 0) {
						y = (N1 - r * N2) / (N1 + r * N2);
						ey = N1 * 4 * r * r * N2 * N2 / ((N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2));
						ey += N2 * 4 * r * r * N1 * N1
								/ ((N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2) * (N1 + r * N2));
						ey = Math.sqrt(ey);
					} else {
						System.out.println("no counts for phi bin " + iAngBin + " " + s);
						y = 0;
						ey = 0;
					}
					// DataFitter fitter;
					x = (iAngBin + 0.5) * 2 * Math.PI / numPhiBins;
					// derivative with respect to N1 is 2*r*N2/(N1+N2)^2
					System.out.println("add point x: "+ x + " y: "+ y + " ex: "+ ex + " ey: "+ ey);
					g.addPoint(x, y, ex, ey);
				}
				DataFitter.fit(f1, g, "Q");
				vals[iKinBin] = f1.parameter(0).value();
				valErrs[iKinBin] = f1.parameter(0).error();

				if ((kinCount[binningType.getBinType()][0][iKinBin]
						+ kinCount[binningType.getBinType()][1][iKinBin]) > 0) {
					double normCount=(kinCount[binningType.getBinType()][0][iKinBin]
							+ kinCount[binningType.getBinType()][1][iKinBin]);
					xVals[iKinBin] = this.meanKin[binningType.binType][iKinBin]
							/ normCount;
					m_normCount[iKinBin]=normCount;
					if(binningType==Binning.RunBinning)
					{
						System.out.print("run bin"+ iKinBin+", meanKin: ");
						System.out.print(this.meanKin[binningType.binType][iKinBin]);
						System.out.println(", normFactor: " + normCount+" xval: "+ xVals[iKinBin]);
					}
					
					
					double wyFactor = this.meanWy[binningType.binType][iKinBin]
							/ normCount;
					this.meanHadMom1[binningType.getBinType()][iKinBin]=this.meanHadMom1[binningType.getBinType()][iKinBin]/normCount;
					this.meanHadMom2[binningType.getBinType()][iKinBin]=this.meanHadMom2[binningType.getBinType()][iKinBin]/normCount;
					this.meanHadTheta1[binningType.getBinType()][iKinBin]=this.meanHadTheta1[binningType.getBinType()][iKinBin]/normCount;
					this.meanHadTheta2[binningType.getBinType()][iKinBin]=this.meanHadTheta2[binningType.getBinType()][iKinBin]/normCount;
					this.meanHadPhi1[binningType.getBinType()][iKinBin]=this.meanHadPhi1[binningType.getBinType()][iKinBin]/normCount;
					this.meanHadPhi2[binningType.getBinType()][iKinBin]=this.meanHadPhi2[binningType.getBinType()][iKinBin]/normCount;
					this.meanQ2[binningType.getBinType()][iKinBin]=this.meanQ2[binningType.getBinType()][iKinBin]/normCount;
					this.meanPt[binningType.getBinType()][iKinBin]=this.meanPt[binningType.getBinType()][iKinBin]/normCount;
					this.meanM[binningType.getBinType()][iKinBin]=this.meanM[binningType.getBinType()][iKinBin]/normCount;
					this.m_pairsWithMatch[binningType.getBinType()][iKinBin]=this.m_pairsWithMatch[binningType.getBinType()][iKinBin]/normCount;
					if(binningType==Binning.RunBinning)
					{
						this.meanHadPairMult[iKinBin]=this.meanHadPairMult[iKinBin]/numEvts[iKinBin];
					}
					
					
					if (doG1P) {
						// should be this A'/C' factor, but couldn't find definition
						wyFactor = 1.0;
					}
					vals[iKinBin] /= wyFactor;
					valErrs[iKinBin] /= wyFactor;
					if(binningType==Binning.RunBinning)
					{
						if(doG1P)
							hAsymDistVsRunG1P.fill(vals[iKinBin]/valErrs[iKinBin]);
						else
							hAsymDistVsRun.fill(vals[iKinBin]/valErrs[iKinBin]);
					}
				
					
				} else {
					System.out.println("no mean vals for " + s);
				}
				System.out.println("iKinbin: " + iKinBin + " val: " + vals[iKinBin] + " err: " + valErrs[iKinBin]
						+ " xval: " + xVals[iKinBin]);
				// should save the graph to make sure it looks ok
				saveGraph(g, s);
			}

			String title = "Asyms_" + baseTitle + binningType.getBinningName();
			if(binningType==Binning.RunBinning)
			{
				for(int i=0;i<xVals.length;i++)
				{
					System.out.println("run binning graph: x " + xVals[i]+ " xErr "+xErrVals[i]+" y: "+vals[i]+ " ey: "+ valErrs[i])
;				}
			}
			
			saveKinGraph(xVals, xErrVals, vals, valErrs, title);
			title ="MeanHadMom1_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanHadMom1[binningType.getBinType()], m_normCount, title, true);
			
			title ="MeanHadMom2_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanHadMom2[binningType.getBinType()], m_normCount, title, true);
			
			title ="MeanHadTheta1_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanHadTheta1[binningType.getBinType()], m_normCount, title, true);
			
			title ="MeanHadTheta2_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanHadTheta2[binningType.getBinType()], m_normCount, title, true);
			
			title ="MeanHadPhi1_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanHadPhi1[binningType.getBinType()], m_normCount, title, true);
			
			title ="MeanHadPhi2_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanHadPhi2[binningType.getBinType()], m_normCount, title,true);
			title ="MeanQ2_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanQ2[binningType.getBinType()], m_normCount, title, true);
			title ="MeanPt_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanPt[binningType.getBinType()], m_normCount, title,true);
			
			title ="MeanM_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, meanM[binningType.getBinType()], m_normCount, title, true);
			
			title ="MeanMatchingFraction_" + baseTitle + binningType.getBinningName();
			saveKinGraph(xVals, xErrVals, m_pairsWithMatch[binningType.getBinType()], m_normCount, title, true);
			
			if(binningType==Binning.RunBinning)
			{
				title ="MeanMult_" + baseTitle + binningType.getBinningName();
				saveKinGraph(xVals, xErrVals, meanHadPairMult, xErrVals, title);
			}
				
			
		}

		// fitHisto.setBinContent(bin, value);
		// fitHisto.setBinError(bin, value);

		////////

	}

	// for the mean values where we estimate the uncertainty from the sample size
	void saveKinGraph(double[] xVals, double[] xValErrs, double[] vals, double[] numCounts, String title,
			boolean roughErrors) {
		// vals can actually be much longer, since the kinetic variables have maxKinBins
		// entries
		double[] valErrs = new double[xVals.length];

		System.out.println("length xval " + xVals.length + " xvalErrs " + xValErrs.length + " vals: " + vals.length
				+ " numCounts: " + numCounts.length + " valErrs: " + valErrs.length);
		for (int i = 0; i < xVals.length; i++) {
			valErrs[i] = Math.sqrt(numCounts[i]) / numCounts[i] * vals[i];
		}

		saveKinGraph(xVals, xValErrs, vals, valErrs, title);
	}

	void saveKinGraph(double[] xVals, double[] xValErrs, double[] vals, double[] valErrs, String title) {
		System.out.println("saving kin graph " + title);

		int length = 0;
		for (int i = 0; i < xVals.length; i++) {

			if (xVals[i] != 0)
				length++;
		}
		System.out.println("skipped " + (xVals.length - length) + " zeros ");
		// to skip zeros...
		double m_xVals[] = new double[length];
		double m_vals[] = new double[length];
		double m_valErrs[] = new double[length];
		double m_xValErrs[] = new double[length];
		int counter = 0;
		for (int i = 0; i < xVals.length; i++) {
			if (xVals[i] != 0) {
				m_xVals[counter] = xVals[i];
				m_xValErrs[counter] = xValErrs[i];
				m_vals[counter] = vals[i];
				m_valErrs[counter] = valErrs[i];
				System.out.println("Setting " + m_xVals[counter] + " [" + counter + "]" + " to " + m_xVals[counter]
						+ ", " + m_xValErrs[counter] + " " + m_vals[counter] + " " + m_valErrs[counter]);
				counter++;
			}
		}

		GraphErrors g = new GraphErrors(title, m_xVals, m_vals, m_xValErrs, m_valErrs);
		g.setTitleY(title);
		g.setTitleX("kinematics");

		EmbeddedCanvas c1 = new EmbeddedCanvas();
		c1.setSize(800, 600);
		// c1.setAxisTitleSize(24);
		// c1.setAxisFontSize(24);
		// PadMargins margins=new PadMargins();
		// margins.setLeftMargin(20);
		// margins.setBottomMargin(20);
		// margins.setRightMargin(20);
		// margins.setTopMargin(20);
		// c1.getPad().setMargins(margins);
		// c1.getPad(0).setMargins(margins);
		String fname = title + ".png";
		c1.update();
		c1.draw(g);
		c1.update();
		c1.save(fname);
		System.out.println("done");

	}

	void saveGraph(GraphErrors g, String title) {
		System.out.println("saving normal graph");
		// Initialize EmbeddedCanvas and divide it
		EmbeddedCanvas c1 = new EmbeddedCanvas();
		c1.setSize(800, 600);
		// c1.setAxisTitleSize(24);
		// c1.setAxisFontSize(24);
		String fname = "fitFor" + title + ".png";
		c1.update();
		c1.draw(g);
		c1.update();
		c1.save(fname);
		System.out.println("done");
	}

	double getWeight(HadronPairData data, double rndSpin) {
		double weight = 1.0;
		double asym = 0.0;
		double sign = 1.0;

		if (rndSpin < 0)
			sign = -1.0;
		weight = 1.0 + sign * asym * Math.sin(data.phiR);
		return weight;

	}

	boolean getHWPInfo(int runNumber) throws Exception {
		// 4239-4326 OUT
		// 4122-4238 IN
		// 3999-4120 OUT
		// 3819-3998 IN
		// 3690-3818 OUT
		// 3479-3551 IN
		// 3243-3475 OUT
		// 3232-3241 IN
		// 3213-3229 OUT
		if (runNumber >= 4239 && runNumber <= 4326)
			return false;
		if (runNumber >= 4122 && runNumber <= 4238)
			return true;
		if (runNumber >= 3999 && runNumber <= 4120)
			return false;
		if (runNumber >= 3819 && runNumber <= 3998)
			return true;
		if (runNumber >= 3690 && runNumber <= 3818)
			return false;
		if (runNumber >= 3479 && runNumber <= 3551)
			return true;
		if (runNumber >= 3243 && runNumber <= 3475)
			return false;
		if (runNumber >= 3232 && runNumber <= 3241)
			return true;
		if (runNumber >= 3213 && runNumber <= 3229)
			return false;
		// uncovered range
		throw new Exception("hwp runnumber " + runNumber + " unknown");

	}

}
