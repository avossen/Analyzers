import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.io.*;
import java.nio.file.*;
import java.io.Serializable;
import java.util.*;

import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.*;

public class PairReader {

	protected H1F phiRResolution;
	protected H1F thetaResolution;
	protected H1F zResolution;
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

	static int numPhiBins = 16;
	// arrays for asymmetry computation. Let's just to pi+pi for now
	// so this is indexed in the kinBin, spin state, phi bin
	protected float[][][] counts;
	protected float[] meanKin;
	protected float[][] kinCount;

	protected ArrayList<Double> phiBins;

	public static void main(String[] args) {

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
		hZ = new H1F("Z", "Z", 100, 0, 1.0);
		hM = new H1F("M", "M", 100, 0, 3.0);
		hXf = new H1F("xF", "xF", 100, -1.0, 1.0);
		hPhiR = new H1F("phiR", "phiR", 100, 0, 2 * Math.PI);
		hPhiH = new H1F("phiH", "phiH", 100, -2 * Math.PI, 2 * Math.PI);
		MVsZResolution = new H2F("MVsZResolution", "MVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		phiRVsZResolution = new H2F("phiRVsZResolution", "phiRVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		hQ2VsX = new H2F("Q2VsX", "Q2VsX", 20, 0.0, 1.0, 20, 0.0, 12);

		counts = new float[Binning.none.numKinBins][2][numPhiBins];
		meanKin = new float[Binning.none.numKinBins];
		// also needed for relative luminosity
		kinCount = new float[Binning.none.numKinBins][2];

		Arrays.fill(counts, (float) 0.0);
		Arrays.fill(meanKin, (float) 0.0);
		Arrays.fill(kinCount, (float) 0.0);

		for (int i = 0; i < numPhiBins; i++) {
			phiBins.add(((i + 1) * 2 * Math.PI / numPhiBins));
		}

	}

	public void savePlots() {
		EmbeddedCanvas can_piPi = new EmbeddedCanvas();
		can_piPi.setSize(1200, 600);
		can_piPi.divide(2, 2);
		can_piPi.setAxisTitleSize(24);
		can_piPi.setAxisFontSize(24);

		can_piPi.setTitleSize(24);
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
		canKin.divide(2, 2);
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
					try {
						// Reading the object from a file
						FileInputStream file = new FileInputStream(filename.toString());
						ObjectInputStream in = new ObjectInputStream(file);

						// Method for deserialization of object
						m_asymData = (AsymData) in.readObject();
						System.out.println("got " + m_asymData.eventData.size() + " hadron pairs");
						for (EventData evtData : m_asymData.eventData) {
							hQ2VsX.fill(evtData.x, evtData.Q2);
							for (HadronPairData pairData : evtData.pairData) {
								hXf.fill(pairData.xF);
								hZ.fill(pairData.z);
								hM.fill(pairData.M);

								if (pairData.z < 0.1 || pairData.xF < 0)
									continue;

								hPhiR.fill(pairData.phiR);
								double weight = 1.0;
								if (pairData.hasMC) {
									weight = getWeight(pairData.matchingMCPair);
								}
								for (Binning binningType : EnumSet.allOf(Binning.class)) {
									int iBin = binningType.getBin(pairData.M, pairData.z, evtData.x);
									int phiBin = binningType.getBin(phiBins, pairData.phiR);
									counts[binningType.binType][evtData.beamPolarization][iBin] += weight;
									kinCount[iBin][evtData.beamPolarization] += weight;
									if (binningType == binningType.MBinning) {
										meanKin[iBin] += pairData.M * weight;

									}
									if (binningType == binningType.ZBinning) {
										meanKin[iBin] += pairData.z * weight;

									}
									if (binningType == binningType.XBinning)
										meanKin[iBin] += evtData.x * weight;
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
		doFits();

		System.out.println("pairs with match " + pairsWithMatch + " pairs without: " + pairsWOMatch + " percentage: "
				+ pairsWithMatch / (float) (pairsWithMatch + pairsWOMatch));
	}

	void doFits() {

		// still no depolarisation factor and beam polarization
		F1D f1 = new F1D("f1", "1+[amp]*sin(x)", -Math.PI, Math.PI);
		f1.setParameter(0, 0.0);

		for (Binning binningType : EnumSet.allOf(Binning.class)) {

			
			
			
			for (int iKinBin = 0; iKinBin < binningType.getNumBins(); iKinBin++) {
				String s = "myFitGraph_" + binningType.getBinningName() + "_bin" + iKinBin;
				GraphErrors g = new GraphErrors(s);

				for (int iAngBin = 0; iAngBin < numPhiBins; iAngBin++) {
					double x = 0;
					double y = 0;
					double ex = 0;
					double ey = 0;

					double r = kinCount[binningType.getBinType()][0] / kinCount[iKin][1];
					double N1 = counts[binningType.getBinType()][0][iAngBin];
					double N2 = counts[binningType.getBinType()][1][iAngBin];
					y = (N1 - r * N2) / (N1 + r * N2);

					// DataFitter fitter;
					x = (iAngBin + 0.5) * 2 * Math.PI / numPhiBins;
					// derivative with respect to N1 is 2*r*N2/(N1+N2)^2
					ey = N1 * 4 * r * r * N2 * N2 / ((N1 + N2) * (N1 + N2) * (N1 + N2) * (N1 + N2));
					ey += N2 * 4 * r * r * N1 * N1 / ((N1 + N2) * (N1 + N2) * (N1 + N2) * (N1 + N2));
					ey = Math.sqrt(ey);

					g.addPoint(x, y, ex, ey);
				}
				DataFitter.fit(f1, g, "Q");
				double amp = f1.parameter(0).value();
				double ampErr=f1.parameter(0).error();
		
				// should save the graph to make sure it looks ok
				saveGraph(g, s);
			}
		}

		// fitHisto.setBinContent(bin, value);
		// fitHisto.setBinError(bin, value);

		////////

	}

	void saveGraph(GraphErrors g, String title) {

		// Initialize EmbeddedCanvas and divide it
		EmbeddedCanvas c1 = new EmbeddedCanvas();
		c1.setSize(1200, 600);
		c1.setAxisTitleSize(24);
		c1.setAxisFontSize(24);
		String fname = "fitFor" + title + ".png";
		c1.draw(g);
		c1.save(fname);
	}

	double getWeight(HadronPairData data) {
		double weight = 1.0;
		return weight;
	}

}
