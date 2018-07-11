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
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;


public class PairReader {

	protected H1F phiRResolution;
	protected H1F thetaResolution;
	protected H1F zResolution;
	protected H1F MResolution;
	
	protected H1F hPhiR;
	protected H1F hPhiH;
	//protected H1F hPhiR;
	
	protected H1F hZ;
	protected H1F hXf;
	
	protected H1F hM;
	
	protected H2F phiRVsZResolution;
	protected H2F MVsZResolution;
	
	
	
	public static void main(String[] args) 
	{
		
		PairReader myReader=new PairReader();
		myReader.initialize();
		myReader.analyze(args);
		myReader.savePlots();
		
		
		
	}
	
	
	public void initialize()
	{
		phiRResolution = new H1F("phiRResolution", "phiRResolution", 100, -1.0, 1.0);
		thetaResolution = new H1F("thetaResolution", "thetaResolution", 100, -1.0, 1.0);
		zResolution = new H1F("zResolution", "zResolution", 100, -1.0, 1.0);
		MResolution = new H1F("MResolution", "MResolution", 100, -1.0, 1.0);
		hZ=new H1F("Z","Z",100,0,1.0);
		hM=new H1F("M","M",100,0,3.0);
		hXf=new H1F("xF","xF",100,-1.0,1.0);
		MVsZResolution = new H2F("MVsZResolution", "MVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
		phiRVsZResolution = new H2F("phiRVsZResolution", "phiRVsZResolution", 20, 0.0, 1.0, 20, -1.0, 1.0);
	}
	public void savePlots()
	{
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
		can_piPi.cd(1);;
		can_piPi.draw(thetaResolution);
		can_piPi.cd(2);;
		can_piPi.draw(zResolution);
		can_piPi.cd(3);;
		can_piPi.draw(MResolution);
		can_piPi.save("resolutions.png");	
		
		EmbeddedCanvas can2D =new EmbeddedCanvas();
		can2D.setSize(1200,600);
		can2D.divide(2,1);
		can2D.cd(0);
		can_piPi.getPad(0).getAxisY().setTitle("phiR Resolution");
		can_piPi.getPad(0).getAxisX().setTitle("z");
		can_piPi.getPad(1).getAxisY().setTitle("M Resolution");
		can_piPi.getPad(0).getAxisX().setTitle("z");
		can2D.draw(this.phiRVsZResolution);
		can2D.cd(1);
		can2D.draw(this.MVsZResolution);
		can2D.save("twoDResolutions.png");
		
		
		EmbeddedCanvas canKin =new EmbeddedCanvas();
		canKin.setSize(1200,600);
		canKin.divide(2, 2);
		canKin.cd(0);
		canKin.getPad(0).getAxisY().setTitle("z");
		canKin.draw(this.hZ);
		canKin.cd(1);
		canKin.getPad(1).getAxisY().setTitle("M");
		canKin.draw(this.hM);
		canKin.getPad(2).getAxisY().setTitle("xF");
		canKin.cd(2);
		canKin.draw(this.hXf);
		canKin.save("diHadKins.png");
		
	}
	
	public void analyze(String[] args)
	{	
		int pairsWithMatch=0;
		int pairsWOMatch=0;
		// TODO Auto-generated method stub
		AsymData m_asymData = null;
		File folder = new File(args[0]);
		File[] listOfFiles = folder.listFiles();
		for (int iF = 0; iF < listOfFiles.length; iF++) {
			if (listOfFiles[iF].isFile()) {
				System.out.println("File " + listOfFiles[iF].getName());
				PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:*.{srn}");

				Path filename = Paths.get(listOfFiles[iF].getName());
				if (matcher.matches(filename)) 
				{
	        // Deserialization
					try
					{   
	            // Reading the object from a file
						FileInputStream file = new FileInputStream(filename.toString());
						ObjectInputStream in = new ObjectInputStream(file);
	             
	            // Method for deserialization of object
						m_asymData = (AsymData)in.readObject();
						System.out.println("got "+ m_asymData.eventData.size() + " hadron pairs");
						for(EventData evtData : m_asymData.eventData)
						{
							for(HadronPairData pairData : evtData.pairData)
							{
								hXf.fill(pairData.xF);
								hZ.fill(pairData.z);
								hM.fill(pairData.M);
								
								if(pairData.z <0.1 || pairData.xF <0)
									continue;
								if(pairData.hasMC)
								{
									pairsWithMatch++;
									phiRResolution.fill((pairData.phiR-pairData.matchingMCPair.phiR)/pairData.phiR);
									
									//System.out.println("pair theta:" + pairData.theta +" matching theta : " + pairData.matchingMCPair.theta );
									thetaResolution.fill((pairData.theta-pairData.matchingMCPair.theta)/pairData.theta);
									zResolution.fill((pairData.z-pairData.matchingMCPair.z)/pairData.z);
									MResolution.fill((pairData.M-pairData.matchingMCPair.M)/pairData.M);
									MVsZResolution.fill(pairData.z, (pairData.M-pairData.matchingMCPair.M)/pairData.M);
									phiRVsZResolution.fill(pairData.z, (pairData.phiR-pairData.matchingMCPair.phiR)/pairData.phiR);
									
									System.out.println("Hadron pair with phiR: "+ pairData.phiR +" has mc partner: "+ pairData.matchingMCPair.phiR);		
								}
								else
								{
									
								pairsWOMatch++;
								}
								//System.out.println("looking at hadron pair with z: "+ pairData.z);			
							}
						}
						in.close();
						file.close();
	             
						System.out.println("Object has been deserialized "); 
	        }	         
	        catch(IOException ex)
	        {
	            System.out.println("IOException is caught");
	        }
	        catch(ClassNotFoundException ex)
	        {
	            System.out.println("ClassNotFoundException is caught");
	        }
	 
	    }
				}
			}
		System.out.println("pairs with match " + pairsWithMatch + " pairs without: " +pairsWOMatch + " percentage: " +pairsWithMatch/(float)(pairsWithMatch+pairsWOMatch));
	}	
	}
