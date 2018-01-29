// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.ArrayList;
import java.util.List;

import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.iterator.MatrixIterator;
import org.la4j.iterator.VectorIterator;
import org.la4j.matrix.DenseMatrix;
import org.la4j.matrix.RowMajorSparseMatrix;
import org.la4j.matrix.sparse.CRSMatrix;
import org.la4j.matrix.sparse.CCSMatrix;
import org.la4j.vector.DenseVector;
import org.la4j.vector.SparseVector;
import org.la4j.vector.dense.BasicVector;

import utils.Utility;

/**
 * References: <br>
 * 1. Abdelghani Bellaachia and Mohammed Al-Dhelaan, Random Walks in Hypergraph, Proceedings of the
 * 2013 International Conference on Applied Mathematics and Computational Methods, pp. 187-194.<br>
 * 2.
 *
 * @author Administrator
 *
 */

public class HyperGraph extends AbstractHyperGraph {
    Matrix H;	// vertex-edge adjacency matrix
    List<HyperEdge> hyperedges;
    Matrix HW;	// weighted vertex-edge adjacency matrix   w(v,e) in the paper


    /**
     * Given vertex-edge weights
     * @param HW
     */
    public HyperGraph(Matrix HW) {
        this(HW, true, false);
    }
    /**
     * Construct a hypergraph given vertex-edge weights and weighting schemes, without given edge weights
     * @param HW
     * @param hyperedgeWeighted
     * @param vertexWeighted
     */
    public HyperGraph(Matrix HW, boolean hyperedgeWeighted, boolean vertexWeighted) {
        this(null, HW, false, hyperedgeWeighted, vertexWeighted);
    }
    /**
     * Given vertex-edge weights and edge weights
     * @param WE
     * @param HW
     */
    public HyperGraph(Vector WE, Matrix HW) {
        this(WE, HW, true, false);
    }
    /**
     * Construct a hypergraph given the vertex-edge weights and edge weights
     * @param WE
     * @param HW
     * @param hyperedgeWeighted
     * @param vertexWeighted
     */
    public HyperGraph(Vector WE, Matrix HW, boolean hyperedgeWeighted, boolean vertexWeighted) {
        this(WE, HW, true, true, false);
    }
    /**
     *
     * @param WE
     * @param HW
     * @param initializeWE
     */
    private HyperGraph(Vector WE, Matrix HW, boolean initializeWE, boolean hyperedgeWeighted, boolean vertexWeighted) {
        assert HW.columns() == WE.length();
        this.HW = HW;
        V = HW.rows(); E = HW.columns();
        this.hyperedgeWeighted = hyperedgeWeighted;
        this.vertexWeighted = vertexWeighted;
        // H
        this.H = HW.copy();
        MatrixIterator itH = null;
        if (HW instanceof DenseMatrix) {
            itH = H.iterator();
        } else {
            if (HW instanceof RowMajorSparseMatrix) {
                itH = ((CRSMatrix) H).nonZeroRowMajorIterator();
            } else {
                itH = ((CCSMatrix) H).nonZeroColumnMajorIterator();
            }
        }
        while (itH.hasNext()) {
            double wv = itH.next();
            if (wv > 0)
                itH.set(1);
        }
        //consolidate hypergraph: calculate hyperedges and HW
        consolidateHG();
        // WE
        if (hyperedgeWeighted) {
            if (initializeWE)	// hyperedge weights are initialized by user
                this.WE = WE;
            else {				// hyperedge weights are calculated according HW
                this.WE = BasicVector.constant(E, 0);
                for (int e = 0; e < E; ++e) {
                    double we = this.getHyperedgeDegree(e, vertexWeighted);
                    this.WE.set(e, we);
                }
            }
        } else {	// not hyperedge-weighted
            this.WE = BasicVector.constant(E, 1);
        }
    }

    /**
     * Given H, WE, WL, hyperedgeWeighted and vertexWeighted either explicitly or in default
     * calculate hyperedges and HW
     */
    private void consolidateHG() {	// 需要修改，考虑weighting
        if (hyperedges != null)
            hyperedges.clear();
        hyperedges = new ArrayList<HyperEdge>(E);
        for (int e = 0; e < E; ++e) {
            HyperEdge he = new HyperEdge(e);
            hyperedges.add(he);
        }
        //
        MatrixIterator it = null;
        if (HW instanceof DenseMatrix) {
            it = HW.iterator();
        } else {
            if (HW instanceof RowMajorSparseMatrix)
                it = ((CRSMatrix) HW).nonZeroRowMajorIterator();
            else
                it = ((CCSMatrix) HW).nonZeroColumnMajorIterator();
        }
        while (it.hasNext()) {
            double w = it.next();
            int v = it.rowIndex();
            int e = it.columnIndex();
            if (w > 0)
                hyperedges.get(e).addVertex(v);
        }
        for (int e = 0; e < E; ++e)
            hyperedges.get(e).consolidate();
    }



    public final Matrix getAdjacencyMatrix(boolean vertexWeighted) {
        if (vertexWeighted)
            return HW;
        else
            return H;
    }



    public final Vector getHyperedgeWeights() {
        return WE;
    }

    /**
     * returns the corresponding value in the hyperedge weight vector WE
     * @param e
     * @return the weight of the e-th hyperedge
     */
    public double getHyperedgeWeight(int e) {
        return WE.get(e);
    }


    /**
     * Equivalent to getAdjacencyMatrix(true).
     * @param vertextWeighted
     * @return
     */
    public final Matrix getVertexWeights() {
        return HW;
    }

    /**
     * Get the weights of the v-th vertext corresponding to all the hyperedges
     * @param v
     * @return a vector of weights of the v-th vertex corresponding to all the hyperedges
     */
    public final Vector getVertexWeight(int v) {
        return HW.getRow(v);
    }

    /**
     * Get the weight of the v-th vertex corresponding to the e-th hyperedge
     * @param v
     * @param e
     * @return the weight of the v-th vertex corresponding to the e-th hyperedge
     */
    public double getVertexWeight(int v, int e) {
        return HW.get(v, e);
    }


    /**
     * vertex/hyperedge degree is different from weight in that degree distinguishes between the weighted and unweighted
     * cases, corresponding to the assignment of <i>vertexWeighted</i>, while weight has only the weighted case
     * @param weighted
     * @return
     */
    @Override
    public final List<Double> getHyperedgeDegrees(boolean vertexWeighted) {
        List<Double> vexEdgeDeg = new ArrayList<Double>(E);
        for (int e = 0; e < E; ++e) {
            vexEdgeDeg.add( getHyperedgeDegree(e, vertexWeighted) );
        }
        return vexEdgeDeg;
    }


    /**
     * Get the degree of the e-th hyperedge, according to the parameter weighted
     * deg(e) = SIGMA_{v incident to e}{w(v,e)}    w(v,e) = HW(v,e) or H(v,e)
     * @param e
     * @param weighted if weighted, w(v,e) = HW(v,e); otherwise w(v,e) = HW(v,e)
     * @return
     */
    @Override
    public double getHyperedgeDegree(int e, boolean vertexWeighted) {
        double deg = 0;
        for (int v : hyperedges.get(e).getVertexIDs()) {
            if (vertexWeighted)
                deg += HW.get(v, e);
            else
                deg++;
        }
        return deg;
    }

    @Override
    public List<Double> getHyperedgeInDegrees(boolean vertexWeighted) {
        // TODO Auto-generated method stub
        return getHyperedgeDegrees(vertexWeighted);
    }
    @Override
    public double getHyperedgeInDegree(int e, boolean vertexWeighted) {
        // TODO Auto-generated method stub
        return getHyperedgeDegree(e, vertexWeighted);
    }

    @Override
    public List<Double> getHyperedgeOutDegrees(boolean vertexWeighted) {
        // TODO Auto-generated method stub
        return getHyperedgeDegrees(vertexWeighted);
    }
    @Override
    public double getHyperedgeOutDegree(int e, boolean vertexWeighted) {
        // TODO Auto-generated method stub
        return getHyperedgeDegree(e, vertexWeighted);
    }

    /**
     * Get the degrees of all the vertex
     * @return
     */
    @Override
    public List<Double> getVertexDegrees(boolean hyperedgeWeighted) {
        List<Double> vecVexDeg = new ArrayList<Double>(V);
        for (int v = 0; v < V; ++V) {
            vecVexDeg.add( getVertexDegree(v, hyperedgeWeighted) );
        }
        return vecVexDeg;
    }

    /**
     *
     * @param v
     * @return
     */
    @Override
    public double getVertexDegree(int v, boolean hyperedgeWeighted) {
        double deg = 0;
        List<Integer> heids = getIncidentHyperedgeIDs(v);
        for (int heid : heids) {
            if (hyperedgeWeighted)
                deg += WE.get(heid);
            else
                deg++; //+= H.get(v, heid);
        }
        return deg;
    }

    @Override
    public List<Double> getTailVertexDegrees(boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        return getVertexDegrees(hyperedgeWeighted);
    }
    @Override
    public double getTailVertexDegree(int v, boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        return getVertexDegree(v, hyperedgeWeighted);
    }
    @Override
    public List<Double> getHeadVertexDegrees(boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        return getVertexDegrees(hyperedgeWeighted);
    }
    @Override
    public double getHeadVertexDegree(int v, boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        return getVertexDegree(v, hyperedgeWeighted);
    }



    public Vector findSinkVertices() {
        Vector sinks = BasicVector.zero(V);
        for (int v = 0; v < V; ++v) {
            if (isSinkVertex(v))
                sinks.set(v, 1);
        }
        return sinks;
    }
    public boolean isSinkVertex(int v) {
        double dv = getVertexDegree(v, hyperedgeWeighted);
        return Utility.isZero(dv);
    }

    public Vector findSinkHyperedges() {
        Vector sinks = BasicVector.zero(E);
        for (int e = 0; e < E; ++e) {
            if (isSinkHyperedge(e))
                sinks.set(e, 1);
        }
        return sinks;
    }
    public boolean isSinkHyperedge(int e) {
        double de = getHyperedgeDegree(e, vertexWeighted);
        return Utility.isZero(de);
    }



    public final List<Integer> getIncidentHyperedgeIDs(int v) {
        List<Integer> heids = new ArrayList<Integer>();
        VectorIterator it = H.iteratorOfRow(v);
        Vector row = H.getRow(v);
        if (row instanceof DenseVector)
            it = row.iterator();
        else
            it = ((SparseVector) row).nonZeroIterator();
        while (it.hasNext()) {
            double w = it.next();
            int he = it.index();
            if (w > 0)
                heids.add(he);
        }
        return heids;
    }

    /**
     * Get the ids of the hyperedges to which the v-th vertex is incident
     * @param v
     * @return
     */
    public final List<HyperEdge> getIncidentHyperedges(int v) {
        List<Integer> heids = getIncidentHyperedgeIDs(v);
        List<HyperEdge> hes = new ArrayList<HyperEdge>(heids.size());
        for (int heid : heids) {
            hes.add(hyperedges.get(heid));
        }
        return hes;
    }

    /**
     * Get the ids of the vertices incident to the e-th hyperedge
     * @param e
     * @return
     */
    public final List<Integer> getIncidentVertices(int e) {
        return hyperedges.get(e).getVertexIDs();
    }

    /**
     * Does not deal with dangling nodes
     * @return
     */
    /**
     * Normalize adjacency matrix HW to transition matrix P
     *
     * 1. neither hyperedge-weighted nor vertex-weighted <br>
     * P = Dv^(-1) * H * De^(-1) * H^(T) <br>
     * where <br>
     * &nbsp;&nbsp;H[u,e] = 1 if e contains u <br>
     * &nbsp;&nbsp;Dv[u] = SIGMA_e{ H[u,e] } = |E(u)| <br>
     * &nbsp;&nbsp;De[e] = SIGMA_v{ H[v,e] } = |V(e)| <br>
     * &nbsp;&nbsp;We = diag(WE) = I <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;WE = constant 1-vector <br>
     * &nbsp;&nbsp;HW = HE <br>
     *
     * 2. hyperedge-weighted but not vertex weighted <br>
     * P = Dv^(-1) * H * We * De^(-1) * H^(T) <br>
     * where <br>
     * &nbsp;&nbsp;H[u,e] = 1 if e contains u <br>
     * &nbsp;&nbsp;Dv[u] = SIGMA_e{ WE[e] * H[u,e] } <br>
     * &nbsp;&nbsp;De[e] = SIGMA_v{ H[v,e] } = |V(e)| <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;HW = HE <br>
     *
     * 3. both hyperedge-weighted and vertex-weighted <br>
     * P = Dv^(-1) * H * We * De^(-1) * HW^(T) <br>
     * where <br>
     * &nbsp;&nbsp;H[u,e] = 1 if e contains u <br>
     * &nbsp;&nbsp;Dv[u] = SIGMA_e{ WE[e] * H[u,e] } <br>
     * &nbsp;&nbsp;De[e] = SIGMA_v{ HW[v,e] }  <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;HW[u,e] > 0, H[u,e] = 1 if e contains u <br>
     *
     * 4. All in one, the following is a general form <br>
     * P = Dv^(-1) * H * We * De^(-1) * HW^(T) <br>
     * where <br>
     * &nbsp;&nbsp;H[u,e] = 1 if e contains u <br>
     * &nbsp;&nbsp;Dv[u] = SIGMA_e{ WE[e] * H[u,e] } <br>
     * &nbsp;&nbsp;De[e] = SIGMA_v{ HW[v,e] } <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;&nbsp;HW[u,e] > 0, H[u,e] = 1 if e contains u <br>
     * &nbsp;If the hypergraph is neither hyperedge-weighted nor vertex-weighted <br>
     * &nbsp;&nbsp;Dv[u] = |E(u)| <br>
     * &nbsp;&nbsp;De[e] = |V(e)| <br>
     * &nbsp;&nbsp;We = I <br>
     * &nbsp;&nbsp;HW = H <br>
     * &nbsp;If the hypergraph is hyperedge-weighted but not vertex-weighted <br>
     * &nbsp;&nbsp;De[e] = |V(e)| <br>
     * &nbsp;&nbsp;HW = H <br>
     * &nbsp;If the hypergraph is both hyperedge-weighted and vertex-weighted <br>
     * @return
     */
    @Override
    public Matrix normalize(boolean hyperedgeWeighted, boolean vertexWeighted) {
        Matrix P = null;
//		int V = this.getNumberOfVertices(), E = this.getNumberOfHyperedges();
        // H and Hw
        Matrix H = this.getAdjacencyMatrix(false);
        Matrix HW = this.getAdjacencyMatrix(vertexWeighted);
        Vector sinkVertices = this.findSinkVertices();
        Vector sinkHyperedges = this.findSinkHyperedges();
//		// Inverse of Dv: d(v) = this.getVertexDegree(v); d(v) = SIGMA_{e}{WE[e]*H[v,e]}
//		Matrix DvI = CRSMatrix.diagonal(V, 1);
//		for (int v = 0; v < V; ++v) {
//			if (Utility.isZero(sinks.get(v))) {	// not sink
//				double dv = this.getVertexDegree(v, hyperedgeWeighted);
//				DvI.set(v, v, 1.0/dv);
//			} else
//				DvI.set(v, v, 0);
//		}
//		// We: diagonal matrix of hyperedge weights; We(e) = WE[e]
//		Matrix We = CRSMatrix.diagonal(E, 1);
//		for (int e = 0; e < E; ++e) {
//			double we = this.getHyperedgeWeight(e);
//			We.set(e, e, we);
//		}
//		// Inverse of De
//		Matrix DeI = CRSMatrix.diagonal(E, 1);
//		for (int e = 0; e < E; ++e) {
//			double de = this.getHyperedgeDegree(e, vertexWeighted);
//			DeI.set(e, e, 1.0/de);
//		}
//		//
//		P = (DvI.multiply(H).multiply(We)).multiply(DeI.multiply(HW.transpose()));
        // reduce the number of multiplications: codes ...
        Vector Dv = BasicVector.zero(V);
        for (int v = 0; v < V; ++v) {
            if (Utility.isZero(sinkVertices.get(v))) {	// not sink vertices
                double dv = this.getVertexDegree(v, hyperedgeWeighted);
                Dv.set(v, dv);
            }
        }
        Vector We = BasicVector.constant(E, 1);
        Vector De = BasicVector.zero(E);
        for (int e = 0; e < E; ++e) {
            if (hyperedgeWeighted) {
                double we = this.getHyperedgeWeight(e);
                We.set(e, we);
            }
            if (Utility.isZero(sinkHyperedges.get(e))) {	// not sink hyperedges
                double de = this.getHyperedgeDegree(e, vertexWeighted);
                De.set(e, de);
            }
        }
        // P = Dv^(-1) * H * We * De^(-1) * HW^(T)
        Matrix P1 = H.copy();
        MatrixIterator itP1 = null;
        if (P1 instanceof DenseMatrix) {
            itP1 = P1.iterator();
        } else {
            if (P1 instanceof RowMajorSparseMatrix) {
                itP1 = ((CRSMatrix) P1).nonZeroRowMajorIterator();
            } else {
                itP1 = ((CCSMatrix) P1).nonZeroColumnMajorIterator();
            }
        }
        // 1. Dv^(-1) * H * We * De^(-1);  P1 = H    Correctness verified
        while (itP1.hasNext()) {
            double wv = itP1.next();
            int v = itP1.rowIndex();
            if (!Utility.isZero(sinkVertices.get(v)))	// sink vertices
                continue;
            wv = wv / Dv.get(v);				// divide each v-th row of P1 by Dv(v)
            int e = itP1.columnIndex();
            if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
                continue;
            wv = wv * We.get(e) / De.get(e);	// multiply each e-th cloumn by We(e) and divide the same column by De(e)
            itP1.set(wv);
        }
        // 2. Dv^(-1) * H * We * De^(-1) * HW;
        P = P1.multiply(HW.transpose());
//		for (int i = 0; i < 1000; ++i) {
//			System.out.printf("sum of row[%d]: %f\n", i, P.getRow(i).fold(Vectors.asSumAccumulator(0)));
//		}
        return P;
    }

    @Override
    public Matrix normalize(boolean H2A, boolean hyperedgeWeighted, boolean vertexWeighted) {
        return normalize(hyperedgeWeighted, vertexWeighted);
    }

    //
//
//	public static void main(String[] args) {
//		// TODO Auto-generated method stub
//		Matrix P = CRSMatrix.unit(4, 4);
//		Matrix P1 = P.copy();
//		Matrix D1 = Basic2DMatrix.diagonal(4, 1);
//		D1.set(0, 0, 1); D1.set(1, 1, 2); D1.set(2, 2, 3); D1.set(3, 3, 4);
//		Matrix D2 = Basic2DMatrix.diagonal(4, 1);
//		D2.set(0, 0, 4); D2.set(1, 1, 4); D2.set(2, 2, 2); D2.set(3, 3, 1);
//		Matrix D3 = Basic2DMatrix.diagonal(4, 1);
//		D3.set(0, 0, 2); D3.set(1, 1, 2); D3.set(2, 2, 4); D3.set(3, 3, 4);
//
//		//
//		MatrixIterator itP1 = null;
//		if (P1 instanceof DenseMatrix) {
//			itP1 = P1.iterator();
//		} else {
//			if (P1 instanceof RowMajorSparseMatrix) {
//				itP1 = ((CRSMatrix) P1).nonZeroRowMajorIterator();
//			} else {
//				itP1 = ((CCSMatrix) P1).nonZeroColumnMajorIterator();
//			}
//		}
//		while (itP1.hasNext()) {
//			double wv = itP1.next();
//			int v = itP1.rowIndex();
//			wv = wv / D1.get(v, v);		    		// divide each v-th row of P1 by Dv(v)
//			int e = itP1.columnIndex();
//			wv = wv * D2.get(e, e) / D3.get(e, e);	// multiply each e-th cloumn by We(e) and divide the same column by De(e)
//			itP1.set(wv);
//		}
//		MatrixUtil.pprint(P1);
//		GaussJordanInverter inverter = new GaussJordanInverter(D1);
//		Matrix D1I = inverter.inverse();
//		inverter = new GaussJordanInverter(D3);
//		Matrix D3I = inverter.inverse();
//		Matrix P2 = D1I.multiply(P.multiply(D2.multiply(D3I)));
////		MatrixUtil.pprint(P2);
//	}
//
    @Override
    public Matrix getOutAdjacencyMatrix(boolean vertexWeighted) {
        // TODO Auto-generated method stub
        return getAdjacencyMatrix(vertexWeighted);
    }

    @Override
    public Matrix getInAdjacencyMatrix(boolean vertexWeighted) {
        // TODO Auto-generated method stub
        return getAdjacencyMatrix(vertexWeighted);
    }
}
