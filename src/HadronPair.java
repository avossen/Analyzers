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

	protected double m_theta;
	protected double phi_h;
	protected double phi_R;

	protected double z;
	protected double xF;

	protected Vector3 Ph;
	protected double pT;
	protected double pTLab;
	protected Vector3 vecQ;
	protected Vector3 vecQLab;
	protected Vector3 vecL;
	protected double m_W;

	protected double m_mass;

	protected MyParticle m_h1;
	protected MyParticle m_h2;
	Vector3 m_breitBoost;
	protected double m_qE;

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
	public HadronPair(MyParticle h1, MyParticle h2, LorentzVector lv_q, LorentzVector lv_l, double W, Vector3 breitBoost) {

		m_h1 = h1;
		m_h2 = h2;
		m_breitBoost = breitBoost;
		m_W = W;
		m_qE = lv_q.e();
		LorentzVector lvQ = new LorentzVector(lv_q);
		LorentzVector lvL = new LorentzVector(lv_l);
		lvQ.boost(breitBoost);
		lvL.boost(breitBoost);
		vecQ = lvQ.vect();
		vecL = lvL.vect();
		vecQLab=new Vector3(lv_q.vect());
		vecQLab.unit();
		compute();

	}

	public void compute() {
		
		LorentzVector pionPair = new LorentzVector(m_h1.px() + m_h2.px(), m_h1.py() + m_h2.py(), m_h1.pz() + m_h2.pz(),
				m_h1.e() + m_h2.e());
		m_mass = pionPair.mass();
		LorentzVector boostedPair = new LorentzVector(pionPair);
		boostedPair.boost(m_breitBoost);
		// we are in the breit frame now, so I guess that the transverse component
		// should already be Pht
		double px = boostedPair.vect().x();
		double py = boostedPair.vect().y();
		
		//I guess the photon is not along the z axis after all
		double phT = Math.sqrt(px * px + py * py);
		pT = phT;
		double otherPt=vecQ.cross(pionPair.vect()).mag();
		//System.out.println("pt " + pT + " or "+otherPT);
		pT=otherPt;
		//mag of cross product should be sin(theta)*|pionPair| (photon vector is unit vector
		pTLab=vecQLab.cross(pionPair.vect()).mag();
		//System.out.println("pt " + pT + " or "+pTLab);
		xF = boostedPair.pz() / m_W;
		z = pionPair.e() / m_qE;

		Vector3 vecQUnit = new Vector3();
		vecQUnit.setMagThetaPhi(1.0, vecQ.theta(), vecQ.phi());
		LorentzVector vh1=new LorentzVector(m_h1.px(),m_h1.py(),m_h1.pz(),m_h1.e());
		LorentzVector vh1T=new LorentzVector(vh1);
		Vector3 pairBoostVect=boostedPair.boostVector();
		pairBoostVect.negative();
		vh1T.boost(pairBoostVect);
		m_theta=Math.acos(vh1T.vect().dot(boostedPair.vect())/(vh1T.vect().mag()*boostedPair.vect().mag()));
		vh1.boost(m_breitBoost);
		LorentzVector vh2=new LorentzVector(m_h2.px(),m_h2.py(),m_h2.pz(),m_h2.e());
		vh2.boost(m_breitBoost);
		
		
		
		
		Vector3 vecRt=new Vector3(vh1.vect());
		vecRt.sub(vh2.vect());
		Vector3 vecPh = boostedPair.vect();
		Vector3 vT = vecQUnit.cross(vecL);
		vT.setMagThetaPhi(1.0, vT.theta(), vT.phi());
		Vector3 vTR = vecQUnit.cross(vecRt);
		Vector3 vTH = vecQUnit.cross(vecPh);
		vTR.setMagThetaPhi(1.0, vTR.theta(), vTR.phi());
		vTH.setMagThetaPhi(1.0, vTH.theta(), vTH.phi());
		double cosPhiR = vT.dot(vTR);
		double cosPhiH = vT.dot(vTH);

		double sinPhiR = vecL.cross(vecRt).dot(vecQUnit);
		double rScale = vecQUnit.cross(vecL).mag() * vecQUnit.cross(vecRt).mag();
		sinPhiR = sinPhiR / rScale;
		double sinPhiH = vecL.cross(vecPh).dot(vecQUnit);
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
		// boostedPair.vect().phi();
	}
}
