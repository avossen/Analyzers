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
import org.jlab.groot.data.GraphErrors;

import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.*;

public class PairReader {

	protected H1F phiRResolution;
	protected H1F thetaResolution;
	protected H1F zResolution;
	protected H1F hY;
	protected H1F hWy;
	protected H1F MResolution;

	protected H1F hPhiR;
	protected H1F hPhiH;
	// protected H1F hPhiR;

	protected H1F hZ;
	protected H1F hXf;

	protected H1F hM;

	protected H2F phiRVsZResolution;
	protected H2F MVsZResolution;
	protected H2F hQ2VsX;

	//is the halfway plate in
	protected boolean hwpIn=false;
	
	static int numPhiBins = 16;
	static int maxKinBins=10;
	// arrays for asymmetry computation. Let's just to pi+pi for now
	// so this is indexed in the kinBin, spin state, phi bin
	protected float[][][][] counts;
	protected float[][][][] countsG1P;
	protected float[][] meanKin;
	protected float[][] meanWy;
	protected float[][][] kinCount;

	protected ArrayList<Double> phiBins;

	public static void main(String[] args) {
		GStyle.getGraphErrorsAttributes().setMarkerStyle(0);
		GStyle.getGraphErrorsAttributes().setMarkerColor(3);
		GStyle.getGraphErrorsAttributes().setMarkerSize(7);
		GStyle.getGraphErrorsAttributes().setLineColor(3);
		GStyle.getGraphErrorsAttributes().setLineWidth(3);
		GStyle.getAxisAttributesX().setTitleFontSize(34);
		GStyle.getAxisAttributesX().setLabelFontSize(10);
		GStyle.getAxisAttributesY().setTitleFontSize(32);
		GStyle.getAxisAttributesY().setLabelFontSize(10);
		
		PairReader myReader = new PairReader();
		myReader.initialize();
		myReader.analyze(args);
		myReader.savePlots();

	}

	public void initialize() {
		phiRResolution = new H1F("phiRResolution", "phiRResolution", 100, -0.3, 0.3);
		thetaResolution = new H1F("thetaResolution", "thetaResolution", 100, -0.3, 0.3);
		zResolution = new H1F("zResolution", "zResolution", 100, -0.3, 0.3);
		MResolution = new H1F("MResolution", "MResolution", 100, -0.3, 0.3);
		hY=new H1F("Y", "Y", 100, -0.3, 1.2);
		hWy=new H1F("Wy", "Wy", 100, -0.3, 1.2);
		hZ = new H1F("Z", "Z", 100, 0, 1.0);
		hM = new H1F("M", "M", 100, 0, 3.0);
		hXf = new H1F("xF", "xF", 100, -1.0, 1.0);
		hPhiR = new H1F("phiR", "phiR", 100, 0, 2 * Math.PI);
		hPhiH = new H1F("phiH", "phiH", 100, -2 * Math.PI, 2 * Math.PI);
		MVsZResolution = new H2F("MVsZResolution", "MVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		phiRVsZResolution = new H2F("phiRVsZResolution", "phiRVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		hQ2VsX = new H2F("Q2VsX", "Q2VsX", 20, 0.0, 1.0, 20, 0.0, 12);
		phiBins=new ArrayList<Double>();
		counts = new float[Binning.numKinBins][2][maxKinBins][numPhiBins];
		countsG1P= new float[Binning.numKinBins][2][maxKinBins][numPhiBins];
		
		meanKin = new float[Binning.numKinBins][maxKinBins];
		meanWy = new float[Binning.numKinBins][maxKinBins];
		// also needed for relative luminosity
		kinCount = new float[Binning.numKinBins][2][maxKinBins];

		//can't do this with multidimensional arrays (but default is already 0)
		//Arrays.fill(counts, (float) 0.0);
		//Arrays.fill(meanKin, (float) 0.0);
		//Arrays.fill(kinCount, (float) 0.0);

		for (int i = 0; i < numPhiBins; i++) {
			System.out.println("adding phi bin "+((i + 1) * 2 * Math.PI / numPhiBins));
			phiBins.add(((i + 1) * 2 * Math.PI / numPhiBins));
		}

	}

	public void savePlots() {
		EmbeddedCanvas can_piPi = new EmbeddedCanvas();
		can_piPi.setSize(1200, 600);
		can_piPi.divide(2, 2);
		//can_piPi.setAxisTitleSize(24);
		//can_piPi.setAxisFontSize(24);

	//	can_piPi.setTitleSize(24);
		can_piPi.cd(0);
		can_piPi.getPad(0).getAxisX().setTitle("phiR Resolution");
		can_piPi.getPad(1).getAxisX().setTitle("theta Resolution");
		can_piPi.getPad(2).getAxisX().setTitle("z Resolution");
		can_piPi.getPad(3).getAxisX().setTitle("M Resolution");
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

		EmbeddedCanvas can2D = new EmbeddedCanvas();
		can2D.setSize(1200, 600);
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
		canAngles.setSize(1200, 600);
		// canAngles.divide(2, 1);
		canAngles.cd(0);
		canAngles.getPad(0).getAxisX().setTitle("phiR");
		canAngles.draw(this.hPhiR);
		// canAngles.getPad(1).getAxisX().setTitle("phiH");
		// canAngles.cd(1);
		// canAngles.draw(this.hPhiH);
		canAngles.save("angles.png");

		EmbeddedCanvas canKin = new EmbeddedCanvas();
		canKin.setSize(1200, 600);
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

	}

	public void analyze(String[] args) {
		int pairsWithMatch = 0;
		int pairsWOMatch = 0;
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
					// Deserialization
					//check run number from filename, so we can determine halfway plate 
					//position
					int runNumber=0;
					StringTokenizer st=new StringTokenizer(listOfFiles[iF].getName(),"_");
				   for(int i=0;st.hasMoreTokens();i++)
				   {
					   String t1=st.nextToken();
					 //  System.out.println("token " + i + ": "+ t1);
					   if(t1.contentEquals("out"))
					   {
						   System.out.println("passed first check");
						   String t2=st.nextToken();
						//   System.out.println("next token " + i + ": "+ t2);
						   if(t2.contentEquals("clas"))
						   {
							  
							   String t3=st.nextToken(".");
							   //need to strip off leading "_" if we use "." token.
							   t3=t3.substring(1, t3.length());
							   //System.out.println("getting run number "+t3);
							   runNumber=Integer.parseInt(t3);
							 
						   }
					   }
				   }
				   //System.out.println("looking at run: " + runNumber);
				   try {
					hwpIn=getHWPInfo(runNumber); 
				   }
				   catch(Exception e)
				   {
					   System.out.println("wrong runnumber");
					   continue;
				   }
					
					try {
						// Reading the object from a file
						ObjectInputStream in=null;
						FileInputStream file=null;
						try 
						{
							file = new FileInputStream(args[0]+"/"+filename.toString());
							in = new ObjectInputStream(file);
						}
						catch(IOException e)
						{
							System.out.println("Caught IO exception: " + e.getMessage());
							throw e;
						}
						// Method for deserialization of object
						m_asymData = (AsymData) in.readObject();
						System.out.println("got " + m_asymData.eventData.size() + " hadron pairs");
						for (EventData evtData : m_asymData.eventData) {
							hQ2VsX.fill(evtData.x, evtData.Q2);
							//kinematic factor W(y) (asymmetry is ~W(y) x e(x)
							double Wy=0;
							if(evtData.y>0)
								Wy=2*evtData.y*Math.sqrt(1-evtData.y);
							else
								continue;
							
							hY.fill(evtData.y);
							hWy.fill(Wy);
							for (HadronPairData pairData : evtData.pairData) {
								hXf.fill(pairData.xF);
								hZ.fill(pairData.z);
								hM.fill(pairData.M);
								//System.out.println("helicity: " + evtData.beamHelicity);
								if (pairData.z < 0.1 || pairData.xF < 0)
									continue;

								hPhiR.fill(pairData.phiR);
								double weight = 1.0;
								if (pairData.hasMC) {
									weight = getWeight(pairData.matchingMCPair);
								}
								int phiBin = Binning.getBin(phiBins, pairData.phiR);
								int phiBinG1P = Binning.getBin(phiBins, pairData.phiH-pairData.phiR);
								//System.out.println("phiR: "+pairData.phiR+ " bin: "+ phiBin);
								for (Binning binningType : EnumSet.allOf(Binning.class)) {
									int iBin = binningType.getBin(pairData.M, pairData.z, evtData.x);
									//int phiBin = binningType.getBin(phiBins, pairData.phiR);
						
								//	System.out.println("kin bin " + binningType.name() + " m: " +pairData.M +" z: "+ pairData.z +" x: " + evtData.x+ " bin: "+ iBin);
				
									if(iBin<0)
									{
										System.out.println("kinematic bin too small " + binningType.name() + " m: " +pairData.M +" z: "+ pairData.z +" x: " + evtData.x);
									}
									if(phiBin<0)
									{
										System.out.println("phi value " + pairData.phiR + " out of bounds");
									}
									
									//do it that way since e.g. for mc data we get -99
									int helicityIndex=0;
									if(evtData.beamHelicity>0)
										helicityIndex=1;
									//flip helicities
									if(hwpIn)
									{
										if(helicityIndex==1)
											helicityIndex=0;
										else
											helicityIndex=1;
									}
									
									
									counts[binningType.binType][helicityIndex][iBin][phiBin] += weight;
									countsG1P[binningType.binType][helicityIndex][iBin][phiBinG1P] += weight*pairData.pTBreit/pairData.M;
									kinCount[binningType.binType][helicityIndex][iBin] += weight;
									meanWy[binningType.binType][iBin] += Wy * weight;
									if (binningType == Binning.MBinning) {
										meanKin[binningType.binType][iBin] += pairData.M * weight;
									
									}
									if (binningType == Binning.ZBinning) {
										meanKin[binningType.binType][iBin] += pairData.z * weight;

									}
									if (binningType == Binning.XBinning)
										meanKin[binningType.binType][iBin] += evtData.x * weight;
								}

								// hPhiH.fill(pairData.ph);
								if (pairData.hasMC) {
									pairsWithMatch++;
									phiRResolution.fill((pairData.phiR - pairData.matchingMCPair.phiR) / pairData.phiR);

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
		doFits(false);
		doFits(true);
		System.out.println("pairs with match " + pairsWithMatch + " pairs without: " + pairsWOMatch + " percentage: "
				+ pairsWithMatch / (float) (pairsWithMatch + pairsWOMatch));
	}

	void doFits(boolean doG1P) {
		float loc_counts[][][][];
		if(doG1P)
			 loc_counts=countsG1P;
		else
			loc_counts=counts;
		// still no depolarisation factor and beam polarization
		//same function for both, e(x) and G1P
		F1D f1 = new F1D("f1", "0+[amp]*sin(x)", 0.0, 2*Math.PI);
		f1.setParameter(0, 0.0);

		for (Binning binningType : EnumSet.allOf(Binning.class)) {

			double vals[]=new double[binningType.getNumBins()];
			double valErrs[]=new double[binningType.getNumBins()];
			double xVals[]=new double[binningType.getNumBins()];
			double xErrVals[]=new double[binningType.getNumBins()];
			String baseTitle;
			if(doG1P)
				baseTitle="G1P_";
			else
				baseTitle="ex_";
			for (int iKinBin = 0; iKinBin < binningType.getNumBins(); iKinBin++) {
				String s = "myFitGraph_"+baseTitle + binningType.getBinningName() + "_bin" + iKinBin;
				GraphErrors g = new GraphErrors(s);

				for (int iAngBin = 0; iAngBin < numPhiBins; iAngBin++) {
					double x = 0;
					double y = 0;
					double ex = 0;
					double ey = 0;

					double r=1.0;
					if(kinCount[binningType.getBinType()][1][iKinBin]>0)
					{
						r = kinCount[binningType.getBinType()][0][iKinBin] / kinCount[binningType.getBinType()][1][iKinBin];
					}
					else
					{
						System.out.println("no counts (r) for phi bin "+iAngBin+" "+s );
					}
					double N1 = loc_counts[binningType.getBinType()][0][iKinBin][iAngBin];
					double N2 = loc_counts[binningType.getBinType()][1][iKinBin][iAngBin];
					if((N1+r*N2)>0)
					{
						y = (N1 - r * N2) / (N1 + r * N2);
						ey = N1 * 4 * r * r * N2 * N2 / ((N1 + N2) * (N1 + N2) * (N1 + N2) * (N1 + N2));
						ey += N2 * 4 * r * r * N1 * N1 / ((N1 + N2) * (N1 + N2) * (N1 + N2) * (N1 + N2));
						ey = Math.sqrt(ey);
					}
					else 
					{
						System.out.println("no counts for phi bin "+iAngBin+" "+s );
						y=0;
						ey=0;
					}
					// DataFitter fitter;
					x = (iAngBin + 0.5) * 2 * Math.PI / numPhiBins;
					// derivative with respect to N1 is 2*r*N2/(N1+N2)^2
					

					g.addPoint(x, y, ex, ey);
				}
				DataFitter.fit(f1, g, "Q");
				vals[iKinBin] = f1.parameter(0).value();
				valErrs[iKinBin]=f1.parameter(0).error();
				if((kinCount[binningType.getBinType()][0][iKinBin]+kinCount[binningType.getBinType()][1][iKinBin])>0)
				{
					xVals[iKinBin]=this.meanKin[binningType.binType][iKinBin]/(kinCount[binningType.getBinType()][0][iKinBin]+kinCount[binningType.getBinType()][1][iKinBin]);
					double wyFactor=this.meanWy[binningType.binType][iKinBin]/(kinCount[binningType.getBinType()][0][iKinBin]+kinCount[binningType.getBinType()][1][iKinBin]);
					if(doG1P)
					{
						//should be this A'/C' factor, but couldn't find definition
						wyFactor=1.0;
					}
					vals[iKinBin]/=wyFactor;
					valErrs[iKinBin]/=wyFactor;
				}
				else
				{
					System.out.println("no mean vals for "+ s);
				}
				// should save the graph to make sure it looks ok
				saveGraph(g, s);
			}
			
			String title="Asyms_"+baseTitle + binningType.getBinningName();
			saveKinGraph(xVals,xErrVals,vals,valErrs, title);
		}

		// fitHisto.setBinContent(bin, value);
		// fitHisto.setBinError(bin, value);

		////////

	}

	void saveKinGraph(double[] xVals, double[] xValErrs, double[] vals,double[] valErrs, String title)
	{
		System.out.println("saving kin graph");
		GraphErrors g=new GraphErrors(title,xVals,vals,xValErrs,valErrs);
		EmbeddedCanvas c1 = new EmbeddedCanvas();
		c1.setSize(1200, 600);
		//c1.setAxisTitleSize(24);
		//c1.setAxisFontSize(24);
		//PadMargins margins=new PadMargins();
		//margins.setLeftMargin(20);
		//margins.setBottomMargin(20);
		//margins.setRightMargin(20);
		//margins.setTopMargin(20);
		//c1.getPad().setMargins(margins);
		//c1.getPad(0).setMargins(margins);
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
		c1.setSize(1200, 600);
		//c1.setAxisTitleSize(24);
		//c1.setAxisFontSize(24);
		String fname = "fitFor" + title + ".png";
		c1.update();
		c1.draw(g);
		c1.update();
		c1.save(fname);
		System.out.println("done");
	}

	double getWeight(HadronPairData data) {
		double weight = 1.0;
		return weight;
	}
	
	boolean getHWPInfo(int runNumber) throws Exception
	{
//		4239-4326     OUT
//		4122-4238     IN
//		3999-4120     OUT
//		3819-3998     IN
//		3690-3818    OUT
//		3479-3551    IN
//		3243-3475   OUT
//		3232-3241   IN
//		3213-3229   OUT
		if(runNumber>=4239 && runNumber <=4326)
			return false;
		if(runNumber>=4122 && runNumber <=4238)
			return true;
		if(runNumber>=3999 && runNumber <=4120)
			return false;
		if(runNumber>=3819 && runNumber <=3998)
			return true;
		if(runNumber>=3690 && runNumber <=3818)
			return false;
		if(runNumber>=3479 && runNumber <=3551)
			return true;
		if(runNumber>=3243 && runNumber <=3475)
			return false;
		if(runNumber>=3232 && runNumber <=3241)
			return true;
		if(runNumber>=3213 && runNumber <=3229)
			return false;
		//uncovered range
		throw new Exception("hwp runnumber "+ runNumber+ " unknown");
		
	}
	

}
