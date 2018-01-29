import hypergraph.DirectedHyperGraph;
import hypergraph.HeterogeneousHyperGraph;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.Vectors;
import org.la4j.vector.dense.BasicVector;
import utils.Utility;

public class HSHMRR {

    private int Np, Nr, Nv, Ns;
    private HeterogeneousHyperGraph prHHG, pvHHG, psHHG;
    private DirectedHyperGraph pDHG, rDHG, vDHG, sDHG;
    private Vector pScore;
    private Matrix pH2A, pA2H, rH2H, vH2H, sH2H, pA2rH, pA2vH, pA2sH, rH2pA, vH2pA, sH2pA;
    private Vector endpH2A, endpA2H, endrH2H, endvH2H, endsH2H, endpA2rH, endpA2vH, endpA2sH, endrH2pA, endvH2pA, endsH2pA;
    private double[][] coefs;

    public HSHMRR(DirectedHyperGraph pDHG, DirectedHyperGraph rDHG, DirectedHyperGraph vDHG, DirectedHyperGraph sDHG,
                  HeterogeneousHyperGraph prHHG, HeterogeneousHyperGraph pvHHG, HeterogeneousHyperGraph psHHG) {
        this(pDHG, rDHG, vDHG, sDHG, prHHG, pvHHG, psHHG,
                BasicVector.constant(pDHG.getNumberOfVertices(), 1.0 / pDHG.getNumberOfVertices()));
    }

    public HSHMRR(DirectedHyperGraph pDHG, DirectedHyperGraph rDHG, DirectedHyperGraph vDHG, DirectedHyperGraph sDHG,
                  HeterogeneousHyperGraph prHHG, HeterogeneousHyperGraph pvHHG, HeterogeneousHyperGraph psHHG,
                  Vector pScore) {
        this(pDHG, rDHG, vDHG, sDHG, prHHG, pvHHG, psHHG,pScore,
                new double[][] {
                        {0.1, 0.38, 0.17, 0.17, 0.17, 0.01},
                        {0.1, 0.9},
                        {0.1, 0.9},
                        {0.1, 0.9}
                });
    }

    public HSHMRR(DirectedHyperGraph pDHG, DirectedHyperGraph rDHG, DirectedHyperGraph vDHG, DirectedHyperGraph sDHG,
                  HeterogeneousHyperGraph prHHG, HeterogeneousHyperGraph pvHHG, HeterogeneousHyperGraph psHHG,
                  Vector pScore, double[][] coefs) {
        this.pDHG = pDHG;
        this.rDHG = rDHG;
        this.vDHG = vDHG;
        this.sDHG = sDHG;
        this.prHHG = prHHG;
        this.pvHHG = pvHHG;
        this.psHHG = psHHG;
        this.pScore = pScore;
        this.coefs = coefs;
        this.Np = pDHG.getNumberOfVertices();
        this.Nr = rDHG.getNumberOfVertices();
        this.Nv = vDHG.getNumberOfVertices();
        this.Ns = sDHG.getNumberOfVertices();
        CalTransitionMatrices();
        ends();
    }

    private void CalTransitionMatrices() {
        System.out.println("######Calculating the transition matrix of paper directed hypergraph######");
        pH2A = pDHG.normalize(true);
        pA2H = pDHG.normalize(false);
        pH2A = pH2A.transpose();
        pA2H = pA2H.transpose();
        System.out.println("######Calculating the transition matrix of author directed hypergraph######");
        rH2H = rDHG.normalize(true);
        rH2H = rH2H.transpose();
        System.out.println("######Calculating the transition matrix of venue directed hypergraph######");
        vH2H = vDHG.normalize(true);
        vH2H = vH2H.transpose();
        System.out.println("######Calculating the transition matrix of institution directed hypergraph######");
        sH2H = sDHG.normalize(true);
        sH2H = sH2H.transpose();
        System.out.println("######Calculating the transition matrix of paper-author heterogeneous hypergraph######");
        Matrix[][] mat = prHHG.normalize();
        pA2rH = mat[0][1];
        rH2pA = mat[1][0];
        pA2rH = pA2rH.transpose();
        rH2pA = rH2pA.transpose();
        System.out.println("######Calculating the transition matrix of paper-venue heterogeneous hypergraph######");
        mat = pvHHG.normalize();
        pA2vH = mat[0][1];
        vH2pA = mat[1][0];
        pA2vH = pA2vH.transpose();
        vH2pA = vH2pA.transpose();
        System.out.println("######Calculating the transition matrix of paper-institution heterogeneous hypergraph######");
        mat = psHHG.normalize();
        pA2sH = mat[0][1];
        sH2pA = mat[1][0];
        pA2sH = pA2sH.transpose();
        sH2pA = sH2pA.transpose();
    }

    private void ends() {
        endpH2A = BasicVector.zero(Np);
        endpA2H = endpH2A.copy();
        endpA2rH = endpH2A.copy();
        endpA2vH = endpH2A.copy();
        endpA2sH = endpH2A.copy();
        for (int i = 0; i < Np; ++i) {
            double sum = pH2A.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endpH2A.set(i, 1);
            sum = pA2H.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endpA2H.set(i, 1);
            sum = pA2rH.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endpA2rH.set(i, 1);
            sum = pA2vH.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endpA2vH.set(i, 1);
            sum = pA2sH.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endpA2sH.set(i, 1);
        }
        endrH2pA = BasicVector.zero(Nr);
        endrH2H = endrH2pA.copy();
        for (int i = 0; i < Nr; ++i) {
            double sum = rH2pA.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endrH2pA.set(i, 1);
            sum = rH2H.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endrH2H.set(i, 1);
        }
        endvH2pA = BasicVector.zero(Nv);
        endvH2H = endvH2pA.copy();
        for (int i = 0; i < Nv; ++i) {
            double sum = vH2pA.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endvH2pA.set(i, 1);
            sum = vH2H.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endvH2H.set(i, 1);
        }
        endsH2pA = BasicVector.zero(Ns);
        endsH2H = endsH2pA.copy();
        for (int i = 0; i < Ns; ++i) {
            double sum = sH2pA.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endsH2pA.set(i, 1);
            sum = sH2H.foldColumn(i, Vectors.asSumAccumulator(0));
            if (Utility.isZero(sum))
                endsH2H.set(i, 1);
        }
    }

    private Vector[] iterate(Vector pA, Vector rH, Vector vH, Vector sH) {
        Vector new_pA, new_pH, new_rH, new_vH, new_sH, tmp;

        // @new_pH <- @pA
        new_pH = pA2H.multiply(pA);
        tmp = BasicVector.constant(Np, 1.0 / Np).multiply(endpA2H.innerProduct(pA));
        new_pH = new_pH.add(tmp);

        // @new_rH <- @rH
        new_rH = rH2H.multiply(rH);
        tmp = BasicVector.constant(Nr, 1.0 / Nr).multiply(endrH2H.innerProduct(rH));
        new_rH = new_rH.add(tmp).multiply(coefs[1][0]);

        // @new_rH <- @pA
        tmp = pA2rH.multiply(pA).multiply(coefs[1][1]);
        new_rH = new_rH.add(tmp);
        tmp = BasicVector.constant(Nr, 1.0 / Nr).multiply(endpA2rH.innerProduct(pA) * coefs[1][1]);
        new_rH = new_rH.add(tmp);

        // @new_vH <- @vH
        new_vH = vH2H.multiply(vH);
        tmp = BasicVector.constant(Nv, 1.0 / Nv).multiply(endvH2H.innerProduct(vH));
        new_vH = new_vH.add(tmp).multiply(coefs[2][0]);

        // @new_vH <- @pA
        tmp = pA2vH.multiply(pA).multiply(coefs[2][1]);
        new_vH = new_vH.add(tmp);
        tmp = BasicVector.constant(Nv, 1.0 / Nv).multiply(endpA2vH.innerProduct(pA) * coefs[2][1]);
        new_vH = new_vH.add(tmp);

        // @new_sH <- @sH
        new_sH = sH2H.multiply(sH);
        tmp = BasicVector.constant(Ns, 1.0 / Ns).multiply(endsH2H.innerProduct(sH));
        new_sH = new_sH.add(tmp).multiply(coefs[3][0]);

        // @new_sH <- pA
        tmp = pA2sH.multiply(pA).multiply(coefs[3][1]);
        new_sH = new_sH.add(tmp);
        tmp = BasicVector.constant(Ns, 1.0 / Ns).multiply(endpA2sH.innerProduct(pA) * coefs[3][1]);
        new_sH = new_sH.add(tmp);

        // @new_pA <- pA
        new_pA = pH2A.multiply(pA);
        tmp = BasicVector.constant(Np, 1.0 / Np).multiply(endpH2A.innerProduct(pA));
        new_pA = new_pA.add(tmp).multiply(coefs[0][0]);

        // @new_pA <- new_pH
        tmp = pH2A.multiply(new_pH).multiply(coefs[0][1]);
        new_pA = new_pA.add(tmp);
        tmp = BasicVector.constant(Np, 1.0 / Np).multiply(endpH2A.innerProduct(new_pH) * coefs[0][1]);
        new_pA = new_pA.add(tmp);

        // @new_pA <- new_rH
        tmp = rH2pA.multiply(new_rH).multiply(coefs[0][2]);
        new_pA = new_pA.add(tmp);
        tmp = BasicVector.constant(Np, 1.0 / Np).multiply(endrH2pA.innerProduct(new_rH) * coefs[0][2]);
        new_pA = new_pA.add(tmp);

        // @new_pA <- new_vH
        tmp = vH2pA.multiply(new_vH).multiply(coefs[0][3]);
        new_pA = new_pA.add(tmp);
        tmp = BasicVector.constant(Np, 1.0 / Np).multiply(endvH2pA.innerProduct(new_vH) * coefs[0][3]);
        new_pA = new_pA.add(tmp);

        // @new_pA <- new_sH
        tmp = sH2pA.multiply(new_sH).multiply(coefs[0][4]);
        new_pA = new_pA.add(tmp);
        tmp = BasicVector.constant(Np, 1.0 / Np).multiply(endsH2pA.innerProduct(new_sH) * coefs[0][4]);
        new_pA = new_pA.add(tmp);

        // @new_pA <- pScore
        new_pA = new_pA.add(pScore.multiply(coefs[0][5]));

        Vector[] result = new Vector[5];
        result[0] = new_pA; result[1] = new_pH; result[2] = new_rH; result[3] = new_vH; result[4] = new_sH;
        return result;
    }

    public Vector[] rank() {
        // initialization
        Vector[] scores = new Vector[5];
        scores[0] = BasicVector.constant(Np, 1.0 / Np);
        scores[1] = scores[0].copy();
        scores[2] = BasicVector.constant(Nr, 1.0 / Nr);
        scores[3] = BasicVector.constant(Nv, 1.0 / Nv);
        scores[4] = BasicVector.constant(Ns, 1.0 / Ns);
        double err = Double.MAX_VALUE;
        int iter = 0;
        while (err >= 1e-5) {
            err = 0;
            Vector[] new_scores = iterate(scores[0], scores[2], scores[3], scores[4]);
            for (int i = 0; i < 5; ++i)
                err += new_scores[i].subtract(scores[i]).fold(Vectors.mkEuclideanNormAccumulator());
            scores = new_scores;
            System.out.println(++iter + "-th iteration, error: " + err);
        }
        return scores;
    }
}
