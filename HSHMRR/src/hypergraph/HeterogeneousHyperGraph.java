// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.iterator.MatrixIterator;
import org.la4j.iterator.VectorIterator;
import org.la4j.matrix.DenseMatrix;
import org.la4j.matrix.RowMajorSparseMatrix;
import org.la4j.matrix.sparse.CCSMatrix;
import org.la4j.matrix.sparse.CRSMatrix;
import org.la4j.vector.DenseVector;
import org.la4j.vector.SparseVector;
import org.la4j.vector.dense.BasicVector;

import utils.Utility;

public class HeterogeneousHyperGraph {
    int V, E;	// number of Vertices and Edges
    int nmode;	// number of entitiy types
    List<Integer> NVs;	// number of vertices for all the modes
    List<Matrix> Hs;	// vertex-edge adjacency matrices for all the modes
    Vector WE;	// weights of hyperedges
    List<HeterogeneousHyperEdge> hyperedges;
    List<Matrix> HWs;	// weighted vertex-edge adjacency matrices for all the modes
    boolean hyperedgeWeighted;
    boolean vertexWeighted;


    /**
     * Construct a heterogeneous hypergraph given the vertex-edge weights and the default weighting scheme
     * i.e. hyperedge weighted while vertex unweighted
     * @param HWs
     */
    public HeterogeneousHyperGraph(List<Matrix> HWs) {
        this(HWs, true, false);
    }
    /**
     * Construct a heterogeneous hypergraph given the vertex-edge weights and the assigned weighting scheme
     * @param HWs
     * @param hyperedgeWeighted
     * @param vertexWeighted
     */
    public HeterogeneousHyperGraph(List<Matrix> HWs, boolean hyperedgeWeighted, boolean vertexWeighted) {
//		this(BasicVector.constant(HWs.get(0).columns(), 1), HWs, hyperedgeWeighted, vertexWeighted);
        //           initializeWE
        this(null, HWs, false, hyperedgeWeighted, vertexWeighted);
    }
    /**
     * Construct a heterogeneous hypergraph given the hyperedge weights and vertex-edge weights,
     * and the default weighting scheme
     * @param WE
     * @param HWs
     */
    public HeterogeneousHyperGraph(Vector WE, List<Matrix> HWs) {
        this(WE, HWs, true, false);
    }
    /**
     * Construct a heterogeneous hypergraph given the hyperedge weights and vertex-edge weights,
     * and the assigned weighting scheme
     * @param WE
     * @param HWs
     * @param hyperedgeWeighted
     * @param vertexWeighted
     */
    public HeterogeneousHyperGraph(Vector WE, List<Matrix> HWs, boolean hyperedgeWeighted, boolean vertexWeighted) {
        //        initializeWE
        this(WE, HWs, true, hyperedgeWeighted, vertexWeighted);
    }

    private HeterogeneousHyperGraph(Vector WE, List<Matrix> HWs, boolean initializeWE, boolean hyperedgeWeighted,
                                    boolean vertexWeighted) {
        this.HWs = HWs;
        nmode = HWs.size();
        NVs = new ArrayList<Integer>(nmode);
        V = 0; E = HWs.get(0).columns();
        Hs = new ArrayList<Matrix>(nmode);
        for (int vmode = 0; vmode < nmode; ++vmode) {
            Matrix HW = HWs.get(vmode);
            int NV = HW.rows(); NVs.add(NV); V += NV;
            Matrix H = HW.copy();
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
//			while (itHW.hasNext()) {
//				itHW.next(); itH.next();
//				int v = itHW.rowIndex();
//				int e = itHW.columnIndex();
//				H.set(v, e, 1);
//			}
            while (itH.hasNext()) {
                double wv = itH.next();
                if (wv > 0)
                    itH.set(1);
            }
            Hs.add(H);
        }
        this.hyperedgeWeighted = hyperedgeWeighted;
        this.vertexWeighted = vertexWeighted;
        // consolidate heterogeneous hypergraph: calculate hyperedges and HW
        consolidateHHG();
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
     * calculate hyperedges
     */
    private void consolidateHHG() {
        // initialize he2vids which in turn is used to generate all the heterogeneous hyperedges
        // heid -- vmode -- vid
        List<List<List<Integer>>> he2vids = new ArrayList<List<List<Integer>>>(E);
        for(int e = 0; e < E; ++e) {
            List<List<Integer>> vids = new ArrayList<List<Integer>>(nmode);
            for (int vmode = 0; vmode < nmode; ++vmode) {
                List<Integer> vidsOfThisMode = new LinkedList<Integer>();
                vids.add(vidsOfThisMode);
            }
            he2vids.add(vids);
        }
        // assign he2vids which in turn is used to generate all the heterogeneous hyperedges
        for (int vmode = 0; vmode < nmode; ++vmode) {
            Matrix HW = HWs.get(vmode);
            MatrixIterator it;
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
                int v = it.rowIndex();		// vertex id
                int e = it.columnIndex();	// hyperedge id
                if (w > 0)
                    he2vids.get(e).get(vmode).add(v);
            }
        }
        // generate all the heterogeneous hyperedges based on he2vids
        if (hyperedges != null)
            hyperedges.clear();
        hyperedges = new ArrayList<HeterogeneousHyperEdge>(E);
        for (int e = 0; e < E; ++e) {
            List<List<Integer>> vids = he2vids.get(e);
            HeterogeneousHyperEdge hhe = new HeterogeneousHyperEdge(e, nmode, vids);
            hyperedges.add(hhe);
        }
    }

//	public void set(Matrix H) {
//		assert H.rows() == V && H.columns() == E;
//		Vector WE = BasicVector.constant(E, 1);
//		set(H, WE);
//	}
//
//	public void set(Matrix H, Vector WE) {
//		assert H.columns() == WE.length() && H.rows() == V && H.columns() == E;
//		Vector WV = BasicVector.constant(V, 1);
//		set(H, WE, WV);
//	}
//
//	public void set(Matrix H, Vector WE, Vector WV) {
//		assert H.rows() == WV.length() && H.columns() == WE.length() && H.rows() == V && H.columns() == E;
//		this.H = H;
//		this.HW = H.copy();
//		this.WE = WE;
//		//
//		consolidateHHG();
//	}

    public final List<Matrix> getAdjacencyMatrices(boolean vertexWeighted) {
        if (vertexWeighted)
            return HWs;
        else
            return Hs;
    }

    public final Matrix getAdjacencyMatrix(int vmode, boolean vertexWeighted) {
        if (vertexWeighted)
            return HWs.get(vmode);
        else
            return Hs.get(vmode);
    }

    public boolean isHyperedgeWeighted() {
        return hyperedgeWeighted;
    }

    public final Vector getHyperedgeWeights() {
        return WE;
    }

    public double getHyperedgeWeight(int e) {
        return WE.get(e);
    }

    public boolean isVertexWeighted() {
        return vertexWeighted;
    }

    public final Matrix getVertexWeights(int vmode) {
        return getAdjacencyMatrix(vmode, true);
    }

    public final Vector getVertexWeights(int v, int vmode) {
        return HWs.get(vmode).getRow(v);
    }

    public double getVertexWeight(int v, int vmode, int e) {
        return HWs.get(vmode).get(v, e);
    }

    public int getNumberOfHyperedges() {
        return E;
    }

    public int getNumberOfVertices() {
        return V;
    }

    public final List<Integer> getIncidentHyperedgeIDs(int v, int vmode) {
        List<Integer> heids = new ArrayList<Integer>();
        Matrix H = Hs.get(vmode);
        VectorIterator it;
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

    public final List<List<Integer>> getIncidentHyperedgeIDs(int v) {
        List<List<Integer>> incident = new ArrayList<List<Integer>>(nmode);
        for (int vmode = 0; vmode < nmode; ++vmode) {
            List<Integer> heids = new ArrayList<Integer>();
            Matrix H = Hs.get(vmode);
            VectorIterator it;
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
            incident.add(heids);
        }
        return incident;
    }

    public final List<HeterogeneousHyperEdge> getIncidentHyperedges(int v, int vmode) {
        List<Integer> heids = getIncidentHyperedgeIDs(v, vmode);
        List<HeterogeneousHyperEdge> hes = new ArrayList<HeterogeneousHyperEdge>(heids.size());
        for (int heid : heids) {
            hes.add(hyperedges.get(heid));
        }
        return hes;
    }

    public final List<List<Integer>> getIncidentVertices(int e) {
        return hyperedges.get(e).getVertexIDs();
    }

    public final List<Integer> getIncidentVertices(int e, int vmode) {
        return hyperedges.get(e).getVertexIDs(vmode);
    }



    public List<Double> getHyperedgeDegrees(int e) {
        return getHyperedgeDegrees(e, vertexWeighted);
    }
    public List<Double> getHyperedgeInDegrees(int e) {
        return getHyperedgeDegrees(e);
    }
    public List<Double> getHyperedgeOutDegrees(int e) {
        return getHyperedgeDegrees(e);
    }

    public List<Double> getHyperedgeDegrees(int e, boolean vertexWeighted) {
        List<Double> degs = new ArrayList<Double>(nmode);
        for (int vmode = 0; vmode < nmode; ++vmode) {
            double deg = getHyperedgeDegree(e, vmode, vertexWeighted);
            degs.add(deg);
        }
        return degs;
    }
    public List<Double> getHyperedgeInDegrees(int e, boolean vertexWeighted) {
        return getHyperedgeDegrees(e, vertexWeighted);
    }
    public List<Double> getHyperedgeOutDegrees(int e, boolean vertexWeighted) {
        return getHyperedgeDegrees(e, vertexWeighted);
    }

    /**
     * Calculate the weighted degree of this hyperedge according to all vertices of different types
     * @param e
     * @return
     */
    public double getHyperedgeDegree(int e) {
        return getHyperedgeDegree(e, vertexWeighted);
    }
    public double getHyperedgeInDegree(int e) {
        return getHyperedgeDegree(e);
    }
    public double getHyperedgeOutDegree(int e) {
        return getHyperedgeDegree(e);
    }

    /**
     * Calculate the weighted/unweighted degree of this hyperedge according to all vertices of different types
     * @param e
     * @param vertexWeighted
     * @return
     */
    public double getHyperedgeDegree(int e, boolean vertexWeighted) {
        double deg = 0;
        for (int vmode = 0; vmode < nmode; ++vmode) {
            deg += getHyperedgeDegree(e, vmode, vertexWeighted);
        }
        return deg;
    }
    public double getHyperedgeInDegree(int e, boolean vertexWeighted) {
        return getHyperedgeDegree(e, vertexWeighted);
    }
    public double getHyperedgeOutDegree(int e, boolean vertexWeighted) {
        return getHyperedgeDegree(e, vertexWeighted);
    }

    /**
     * Calculate the weighted degree of this hyperedge according to the vertices of a specific type (i.e. vmode)
     * @param e
     * @param vmode
     * @return
     */
    public double getHyperedgeDegree(int e, int vmode) {
        return getHyperedgeDegree(e, vmode, vertexWeighted);
    }
    public double getHyperedgeInDegree(int e, int vmode) {
        return getHyperedgeDegree(e, vmode);
    }
    public double getHyperedgeOutDegree(int e, int vmode) {
        return getHyperedgeDegree(e, vmode);
    }

    /**
     * Calculate the weighted/unweighted degree of this hyperedge according to the vertices of a specific type (i.e. vmode)
     * @param e
     * @param vmode
     * @param vertexWeighted
     * @return
     */
    public double getHyperedgeDegree(int e, int vmode, boolean vertexWeighted) {
        double deg = 0;
        for (int v : hyperedges.get(e).getVertexIDs(vmode)) {
            if (vertexWeighted)
                deg += HWs.get(vmode).get(v, e);
            else
                deg++;
        }
        return deg;
    }
    public double getHyperedgeInDegree(int e, int vmode, boolean vertexWeighted) {
        return getHyperedgeDegree(e, vmode, vertexWeighted);
    }
    public double getHyperedgeOutDegree(int e, int vmode, boolean vertexWeighted) {
        return getHyperedgeDegree(e, vmode, vertexWeighted);
    }



    public double getVertexDegree(int v, int vmode) {
        return getVertexDegree(v, vmode, hyperedgeWeighted);
    }
    public double getTailVertexDegree(int v, int vmode) {
        return getVertexDegree(v, vmode);
    }
    public double getHeadVertexDegree(int v, int vmode) {
        return getVertexDegree(v, vmode);
    }

    public double getVertexDegree(int v, int vmode, boolean hyperedgeWeighted) {
        double deg = 0;
        List<Integer> heids = getIncidentHyperedgeIDs(v, vmode);
        for (int heid : heids) {
            if (hyperedgeWeighted)
                deg += WE.get(heid);
            else
                deg++; //+= H.get(v, heid);
        }
        return deg;
    }
    public double getTailVertexDegree(int v, int vmode, boolean hyperedgeWeighted) {
        return getVertexDegree(v, vmode, hyperedgeWeighted);
    }
    public double getHeadVertexDegree(int v, int vmode, boolean hyperedgeWeighted) {
        return getVertexDegree(v, vmode, hyperedgeWeighted);
    }

    /**
     * Does not deal with dangling nodes
     * @return
     */
    public Matrix[][] normalize() {
        return normalize(hyperedgeWeighted, vertexWeighted);
    }
    /**
     * Normalize all the adjacency matrices HWs to transition matrices Ps <br>
     * For each adjacency matrix HW, the corresponding transition matrix P is calculated as follows: <br>
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
     * P = Dv^(-1) * H * We * De^(-1) * HW <br>
     * where <br>
     * &nbsp;&nbsp;H[u,e] = 1 if e contains u <br>
     * &nbsp;&nbsp;Dv[u] = SIGMA_e{ WE[e] * H[u,e] } <br>
     * &nbsp;&nbsp;De[e] = SIGMA_v{ HW[v,e] }  <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;HW[u,e] > 0, H[u,e] = 1 if e contains u <br>
     *
     * 4. All in one, the following is a general form <br>
     * P = Dv^(-1) * H * We * De^(-1) * HW <br>
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
     *
     * @return NormHWs is a set of transition matrices between each two types of entities. <br>
     * Each NormHW[j][i] is the matrix of transition probabilities of jumping from a j-th type entity to an i-th type entity.
     */
    public Matrix[][] normalize(boolean hyperedgeWeighted, boolean vertexWeighted) {
        Matrix[][] NormHWs = new Matrix[nmode][nmode];
        //
        // Each NormHW[j][i] is the matrix of transition probabilities of jumping from
        // a j-th type entity to an i-th type entity.
        for (int j = 0; j < nmode; ++j) {
            Vector sinkVertices = this.findSinkVertices(j);
            for (int i = 0; i < nmode; ++i) {
                if (i == j) continue;
                System.out.printf("Normalizing Heterogeneous Hypergraph NormHWs[%d][%d] ...\n", j, i);
                Vector sinkHyperedges = findSinkHyperedges(i);
                // NormHWs[j][i]: how the j-th mode influences the i-th mode
                Matrix P = null;
                // H and Hw
                Matrix H = this.getAdjacencyMatrix(j, false);
                Matrix HW = this.getAdjacencyMatrix(i, vertexWeighted);
//				// Inverse of Dv
//				int NV = NVs.get(j);
//				Matrix DvI = CRSMatrix.diagonal(NV, 1);
//				for (int v = 0; v < NV; ++v) {
//					if (Utility.isZero(sinkVertices.get(v))) {	// not sink
//						double dv = this.getVertexDegree(v, j, hyperedgeWeighted);
//						DvI.set(v, v, 1.0 / dv);
//					} else
//						DvI.set(v, v, 0);
//				}
//				// We: diagonal matrix of hyperedge weights
//				Matrix We = CRSMatrix.diagonal(E, 1);
//				for (int e = 0; e < E; ++e) {
//					if (hyperedgeWeighted) {
//						double we = this.getHyperedgeWeight(e);
//						We.set(e, e, we);
//					}
//				}
//				// Inverse of De
//				Matrix DeI = CRSMatrix.diagonal(E, 1);
//				for (int e = 0; e < E; ++e) {
//					if (Utility.isZero(sinkHyperedges.get(e))) {	// not sink
//						double de = this.getHyperedgeDegree(e, i, vertexWeighted);
//						DeI.set(e, e, 1.0 / de);
//					} else
//						DeI.set(e, e, 0);
//				}
//				//
//				P = (DvI.multiply(H).multiply(We)).multiply(DeI.multiply(HW.transpose()));
                // reduce the number of multiplications: codes ...
                int NV = NVs.get(j);
                Vector Dv = BasicVector.zero(NV);
                for (int v = 0; v < NV; ++v) {
                    if (Utility.isZero(sinkVertices.get(v))) {	// not sink vertices
                        double dv = this.getVertexDegree(v, j, hyperedgeWeighted);
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
                        double de = this.getHyperedgeDegree(e, i, vertexWeighted);
                        De.set(e, de);
                    }
                }
                // P = Dv^(-1) * H * We * De^(-1) * HW^(T)
                // 1. Dv^(-1) * H * We * De^(-1);  P1 = H    Correctness verified
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
                while (itP1.hasNext()) {
                    double wv = itP1.next();
                    int v = itP1.rowIndex();
                    if (!Utility.isZero(sinkVertices.get(v)))	// sink vertices
                        continue;
                    wv = wv / Dv.get(v);				// divide each v-th row of P1 by Dv(v)
                    int e = itP1.columnIndex();
                    if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
                        continue;
                    wv = wv * We.get(e) / De.get(e);	// multiply each e-th cloumn by We(e)
                    itP1.set(wv);
                }
                // 2. Dv^(-1) * H * We * De^(-1) * HW;
//				Matrix P2 = HW.copy();
//				MatrixIterator itP2 = null;
//				if (P2 instanceof DenseMatrix) {
//					itP2 = P2.iterator();
//				} else {
//					if (P2 instanceof RowMajorSparseMatrix) {
//						itP2 = ((CRSMatrix) P1).nonZeroRowMajorIterator();
//					} else {
//						itP2 = ((CCSMatrix) P1).nonZeroColumnMajorIterator();
//					}
//				}
//				while (itP2.hasNext()) {
//					double wv = itP2.next();
//					int e = itP1.columnIndex();
//					if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
//						continue;
//					wv = wv / De.get(e);	// divide each e-th column by De(e)
//					itP2.set(wv);
//				}
                P = P1.multiply(HW.transpose());
                //
//				for (int k = 0; k < 1000; ++k) {
//					System.out.printf("sum of row[%d]: %f\n", k, P.getRow(k).fold(Vectors.asSumAccumulator(0)));
//				}
                NormHWs[j][i] = P;
            }
        }
        return NormHWs;
    }



    public List<Vector> findSinkVertices() {
        List<Vector> sinks = new ArrayList<Vector>(nmode);
        for (int vmode = 0; vmode < nmode; ++vmode) {
            sinks.add( findSinkVertices(vmode) );
        }
        return sinks;
    }
    public Vector findSinkVertices(int vmode) {
        int NV = NVs.get(vmode);
        Vector sinks = BasicVector.zero(NV);
        for (int v = 0; v < NV; ++v) {
            if (this.isSinkVertex(v, vmode))
                sinks.set(v, 1);
        }
        return sinks;
    }
    public boolean isSinkVertex(int v, int vmode) {
        double dv = this.getVertexDegree(v, vmode, hyperedgeWeighted);
        return Utility.isZero(dv);
    }

    public List<Vector> findSinkHyperedges() {
        List<Vector> sinks = new ArrayList<Vector>(nmode);
        for (int vmode = 0; vmode < nmode; ++vmode) {
            sinks.add( findSinkHyperedges(vmode) );
        }
        return sinks;
    }
    public Vector findSinkHyperedges(int vmode) {
        Vector sinks = BasicVector.zero(E);
        for (int e = 0; e < E; ++e) {
            if (this.isSinkHyperedge(e, vmode))
                sinks.set(e, 1);
        }
        return sinks;
    }
    public boolean isSinkHyperedge(int e, int vmode) {
        double de = this.getHyperedgeDegree(e, vmode, vertexWeighted);
        return Utility.isZero(de);
    }
}
