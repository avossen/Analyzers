import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;
import NovelFitters.MyParticle;

public class HadronPair {
	boolean hasMatchingMC;
	HadronPair matchingMC;
	HadronPairData matchingMCData;
	protected double m_theta;
	protected double phi_h;
	protected double phi_R;

	protected double z;
	protected double xF;
	protected double missingMass;

	protected Vector3 Ph;
	protected double pT;
	protected double pTLab;
	protected Vector3 vecQ;
	protected Vector3 vecQLab;
	protected LorentzVector lvQLab;
	protected Vector3 vecL;
	protected double m_W;

	protected double m_mass;

	protected MyParticle m_h1;
	protected MyParticle m_h2;
	Vector3 m_breitBoost;
	protected double m_qE;

	public double getMissingMass() {
		return missingMass;
	}

	public double getTheta() {
		return m_theta;
	}

	public double getPhiH() {
		return phi_h;
	}

	public double getPhiR() {
		return phi_R;
	}

	public double getZ() {
		return z;
	}

	public double getXf() {
		return xF;
	}

	public Vector3 geetPh() {
		return Ph;
	}

	public double getPt() {
		return pT;
	}

	public double getPtLab() {
		return pTLab;
	}

	public double getMass() {
		return m_mass;
	}

	// takes virtual photon and Was parameter
	public HadronPair(MyParticle h1, MyParticle h2, LorentzVector lv_q, LorentzVector lv_l, double W,
			Vector3 breitBoost) {
		hasMatchingMC=false;	
		m_h1 = h1;
		m_h2 = h2;
		m_breitBoost = new Vector3(breitBoost);
		m_W = W;
		m_qE = lv_q.e();
		LorentzVector lvQ = new LorentzVector(lv_q);
		lvQLab = new LorentzVector(lv_q);
		LorentzVector lvL = new LorentzVector(lv_l);
		lvQ.boost(breitBoost);
		lvL.boost(breitBoost);
		vecQ = lvQ.vect();
		vecL = lvL.vect();
		vecQLab = new Vector3(lv_q.vect());
		vecQLab.unit();
		compute();

	}

	public void compute() {

		LorentzVector pionPair = new LorentzVector(m_h1.px() + m_h2.px(), m_h1.py() + m_h2.py(), m_h1.pz() + m_h2.pz(),m_h1.e() + m_h2.e());
		double m_p = 0.938;
		LorentzVector protonLab = new LorentzVector();
		protonLab.setPxPyPzM(0.0, 0.0, 0.0, m_p);
		LorentzVector tmp = new LorentzVector(protonLab);
	//	tmp.add(protonLab);
		tmp.add(lvQLab);
		tmp.sub(pionPair);
		missingMass = tmp.mass();

		m_mass = pionPair.mass();
		LorentzVector boostedPair = new LorentzVector(pionPair);
		boostedPair.boost(m_breitBoost);
		Vector3 vecQUnit = new Vector3();
		vecQUnit.setMagThetaPhi(1.0, vecQ.theta(), vecQ.phi());
		Vector3 vecQLabUnit = new Vector3();
		vecQLabUnit.setMagThetaPhi(1.0, vecQLab.theta(), vecQLab.phi());
		//double otherPt = vecQ.cross(pionPair.vect()).mag();
		double otherPt = vecQUnit.cross(boostedPair.vect()).mag();
		// System.out.println("pt " + pT + " or "+otherPT);
		pT = otherPt;
		// mag of cross product should be sin(theta)*|pionPair| (photon vector is unit
		// vector
		pTLab = vecQLabUnit.cross(pionPair.vect()).mag();
		//System.out.println("pt " + pT + " or "+pTLab);
		xF = boostedPair.pz() / m_W;
		z = pionPair.e() / m_qE;
		//should be in lab system, so that the proton momentum 
		//doesn't have to be accounted for
		double z1=m_h1.e()/m_qE;
		double z2=m_h2.e()/m_qE;
	
		 //z1=1.0;
		// z2=1.0;
		
		
		LorentzVector vh1 = new LorentzVector(m_h1.px(), m_h1.py(), m_h1.pz(), m_h1.e());
		LorentzVector vh1T = new LorentzVector(vh1);
		Vector3 pairBoostVect = boostedPair.boostVector();
		pairBoostVect.negative();
		vh1T.boost(pairBoostVect);
		m_theta = Math.acos(vh1T.vect().dot(boostedPair.vect()) / (vh1T.vect().mag() * boostedPair.vect().mag()));
		vh1.boost(m_breitBoost);
		LorentzVector vh2 = new LorentzVector(m_h2.px(), m_h2.py(), m_h2.pz(), m_h2.e());
		vh2.boost(m_breitBoost);

		//so vecR is now in the Breit frame
		Vector3 vecR = new Vector3(vh1.vect());
		//scale by 1/z1
		vecR.setMagThetaPhi(vh1.vect().mag()/z1, vh1.vect().theta(),vh2.vect().phi());
		Vector3 scaledVh2=new Vector3();
		scaledVh2.setMagThetaPhi(vh2.vect().mag()/z2, vh2.vect().theta(), vh2.vect().theta());
		vecR.sub(scaledVh2);
		//get vecRT:
		Vector3 vecRt=new Vector3();
		Vector3 RAlongQ=new Vector3();
		
		RAlongQ.setMagThetaPhi(vecR.dot(vecQUnit), vecQUnit.theta(), vecQUnit.phi());
		vecRt=vecR;
		//The transverse part is the original vector minus the one that 
		//is along q
	//	System.out.println("lenght of r along q : " + RAlongQ.mag());
		vecRt.sub(RAlongQ);
	//	System.out.println("length of R " + vecRt.mag() + " phi " + vecRt.phi() + " theta: " + vecRt.theta());
		Vector3 vecPh = new Vector3(boostedPair.vect());
		
		Vector3 PtAlongQ=new Vector3();
		PtAlongQ.setMagThetaPhi(vecR.dot(vecQUnit), vecQUnit.theta(), vecQUnit.phi());
		Vector3 vecPhT=new Vector3(vecPh);
		vecPhT.sub(PtAlongQ);
		
		Vector3 vT = new Vector3(vecQUnit.cross(vecL));
		//vT.setMagThetaPhi(1.0, vT.theta(), vT.phi());
		vT.unit();
		Vector3 vTR = new Vector3(vecQUnit.cross(vecRt));
		Vector3 vTH = new Vector3(vecQUnit.cross(vecPhT));
		vTR.unit();
		vTH.unit();
		double cosPhiR = vT.dot(vTR);
		double cosPhiH = vT.dot(vTH);

		double sinPhiR = vecL.cross(vecRt).dot(vecQUnit);
		///
		//the scaling for cosPhiR is not necessary anymore since 
		//for that quantity we already operate 
		//with unit vectors above, scaling again would lead to the wrong value
		//
		double rScale = vecQUnit.cross(vecL).mag() * vecQUnit.cross(vecRt).mag();
		sinPhiR = sinPhiR / rScale;
		double sinPhiH = vecL.cross(vecPhT).dot(vecQUnit);
		double hScale = vecQUnit.cross(vecL).mag() * vecQUnit.cross(vecPh).mag();
		sinPhiH = sinPhiH / hScale;

		phi_h = Math.acos(cosPhiH);
		if (sinPhiH < 0.0) {
			phi_h = 2 * Math.PI - phi_h;
		}
		phi_R = Math.acos(cosPhiR);
		if (sinPhiR < 0.0) {
			phi_R = 2 * Math.PI - phi_R;
			
		}
		//new way: 
		// turns out it is the same as the old way with the exception that
		//the angle is in the range -pi - pi instead of 0 to 2pi 
		//
		
		double phi_R2Sign=(vecQUnit.cross(vecL)).dot(vecRt)/Math.abs((vecQUnit.cross(vecL)).dot(vecRt));
		double phi_R2=(vecQUnit.cross(vecL)).dot(vecQUnit.cross(vecRt));
		phi_R2=Math.acos(phi_R2/((vecQUnit.cross(vecL)).mag()*(vecQUnit.cross(vecRt)).mag()));
		phi_R2=phi_R2*phi_R2Sign;
		
		
		//System.out.println("Old phiR: " + phi_R+ " new phiR: " + phi_R2);
	//	System.out.println("old sign " + sinPhiR + " new sign: "+ (vecQUnit.cross(vecL)).dot(vecRt)/Math.abs((vecQUnit.cross(vecL)).dot(vecRt)));
		// boostedPair.vect().phi();
	}
}
