// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.javatuples.Pair;
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

import utils.RealAccumulator;
import utils.Utility;

/**
 An examplemplar paper citation hypergraph is illustrated below.<br>
 edge 1: p1 --> p2, p3, p4 <br>
 edge 2: p2 --> p4 <br>
 edge 3: p3 --> p2, p4 <br>
 E = 3, V = 4 <br>
 -           e1  e2  e3                      e1  e2  e3 <br>
 Htail = p1  1   0   0           Hhead = p1  0   0   0  <br>
 p2  0   1   0                   p2  1   0   1  <br>
 p3  0   0   1                   p3  1   0   0  <br>
 p4  0   0   0                   p4  1   1   1  <br>
 For a directed hypergraph, a vertex has two degree, out degree as the weighted sum of the hyperedges to which
 is is incident as tail vertices, and in degree as the weighted sum of the hyperedges to which it is incident as
 head vertices. <br>
 So one one hand, a vertex v might not be head vertices for any hyperedge, then indeg(v) = 0, e.g. vertex p1 in
 the above example, which means p1 is not cited by any paper; on the other hand, a vertex v might not be tail
 vertices for any hyperedge, the outdeg(v) = 0, e.g. p4 in the above example,  <br>
 Considering hyperedges, a hyperedge may have either empty tail or head vertex set. If a hyperedge has either an
 empty tail or head vertex set, there are two ways to deal with the situation, which actually have the same result.
 (M1) We can set the corresponding head or tail vertex set to empty too and thus remove this hyperedge without
 incurring any computation mistakes. In this way, we allow only zero degrees for tail/head vertices and all hyperedges
 must have non-zero in/out degrees, just as the above example shows. (M2) We can also keep these hyperedges which
 have empty tail or head vertices. Thus in this way, either vertex or hyperedge may have zero degrees, i.e. zero in/out
 degrees for hyperedges and zero degrees for tail/head vertices, just as the equivalent example below shows.<br>
 -           e1  e2  e3  e4                      e1  e2  e3  e4 <br>
 Htail = p1  1   0   0   0           Hhead = p1  0   0   0   0  <br>
 p2  0   1   0   0                   p2  1   0   1   0  <br>
 p3  0   0   1   0                   p3  1   0   0   0  <br>
 p4  0   0   0   1                   p4  1   1   1   0  <br>
 * @author Administrator
 *
 */

public class DirectedHyperGraph extends AbstractHyperGraph {

    Matrix Htail;	// vertex-edge adjacency matrix
    Matrix Hhead;	// vertex-edge adjacency matrix
    List<DirectedHyperEdge> hyperedges;
    Matrix HWtail;	// weighted vertex-edge adjacency matrix   w(v,e) in the paper
    Matrix HWhead;	// weighted vertex-edge adjacency matrix   w(v,e) in the paper
    boolean allowSelfLoop;

    /**
     * Given vertex-edge weights
     * @param HWtail
     * @param HWhead
     */
    public DirectedHyperGraph(Matrix HWtail, Matrix HWhead, boolean allowSelfLoop) {
        this(HWtail, HWhead, allowSelfLoop, true, false);	// by default: hyperedgeWeighted = true; vertexWeighted = false
    }
    public DirectedHyperGraph(Matrix HWtail, Matrix HWhead) {
        this(HWtail, HWhead, true);	 // by default: allowSelfLoop = true;
    }

    /**
     * Construct a hypergraph given vertex-edge weights and weighting schemes, without given edge weights
     * @param HWtail
     * @param HWhead
     * @param hyperedgeWeighted
     * @param vertexWeighted
     */
    public DirectedHyperGraph(Matrix HWtail, Matrix HWhead, boolean allowSelfLoop, boolean hyperedgeWeighted, boolean vertexWeighted) {
        this(null, HWtail, HWhead, allowSelfLoop, false, hyperedgeWeighted, vertexWeighted);
    }
    public DirectedHyperGraph(Matrix HWtail, Matrix HWhead, boolean hyperedgeWeighted, boolean vertexWeighted) {
        this(HWtail, HWhead, true, hyperedgeWeighted, vertexWeighted);	//by default: allowSelfLoop = true
    }
    /**
     * Given vertex-edge weights and edge weights
     * @param WE
     * @param HWtail
     * @param HWhead
     */
    public DirectedHyperGraph(Vector WE, Matrix HWtail, Matrix HWhead, boolean allowSelfLoop) {
        this(WE, HWtail, HWhead, allowSelfLoop, true, false);  // by default: hyperedgeWeighted = true; vertexWeighted = false
    }
    public DirectedHyperGraph(Vector WE, Matrix HWtail, Matrix HWhead) {
        this(WE, HWtail, HWhead, true);  // by default: allowSelfLoop = true;
    }
    /**
     * Construct a hypergraph given the vertex-edge weights and edge weights
     * @param WE
     * @param HWtail
     * @param HWhead
     * @param hyperedgeWeighted
     * @param vertexWeighted
     */
    public DirectedHyperGraph(Vector WE, Matrix HWtail, Matrix HWhead, boolean allowSelfLoop, boolean hyperedgeWeighted, boolean vertexWeighted) {
        this(WE, HWtail, HWhead, allowSelfLoop, true, hyperedgeWeighted, vertexWeighted);  // by default: initializeWE = true
    }
    public DirectedHyperGraph(Vector WE, Matrix HWtail, Matrix HWhead, boolean hyperedgeWeighted, boolean vertexWeighted) {
        this(WE, HWtail, HWhead, true, hyperedgeWeighted, vertexWeighted);  // by default: allowSelfLoop = true;
    }
    /**
     *
     * @param WE
     * @param HWtail
     * @param HWhead
     * @param initializeWE
     */
    private DirectedHyperGraph(Vector WE, Matrix HWtail, Matrix HWhead, boolean allowSelfLoop, boolean initializeWE,
                               boolean hyperedgeWeighted, boolean vertexWeighted) {
        this.allowSelfLoop = allowSelfLoop;
//		assert HWtail.columns() == WE.length();
        this.HWtail = HWtail; this.HWhead = HWhead;
        V = HWtail.rows(); E = HWtail.columns();
        this.hyperedgeWeighted = hyperedgeWeighted;
        this.vertexWeighted = vertexWeighted;
        // Htail
        this.Htail = HWtail.copy();
        MatrixIterator itHtail = null;
        if (Htail instanceof DenseMatrix) {
            itHtail = Htail.iterator();
        } else {
            if (Htail instanceof RowMajorSparseMatrix) {
                itHtail = ((CRSMatrix) Htail).nonZeroRowMajorIterator();
            } else {
                itHtail = ((CCSMatrix) Htail).nonZeroColumnMajorIterator();
            }
        }
//		while (itHWtail.hasNext()) {
//			itHWtail.next(); itHtail.next();
//			int v = itHWtail.rowIndex();
//			int e = itHWtail.columnIndex();
//			Htail.set(v, e, 1);
//		}
        while (itHtail.hasNext()) {
            double w = itHtail.next();
            if (w > 0)
                itHtail.set(1);
        }
        // Hhead
        this.Hhead = HWhead.copy();
        MatrixIterator itHhead = null;
        if (Hhead instanceof DenseMatrix) {
            itHhead = Hhead.iterator();
        } else {
            if (Hhead instanceof RowMajorSparseMatrix) {
                itHhead = ((CRSMatrix) Hhead).nonZeroRowMajorIterator();
            } else {
                itHhead = ((CCSMatrix) Hhead).nonZeroColumnMajorIterator();
            }
        }
//		while (itHWhead.hasNext()) {
//			itHWhead.next(); itHhead.next();
//			int v = itHWhead.rowIndex();
//			int e = itHWhead.columnIndex();
//			Hhead.set(v, e, 1);
//		}
        while (itHhead.hasNext()) {
            double wv = itHhead.next();
            if (wv > 0)
                itHhead.set(1);
        }
        //consolidate hypergraph: calculate hyperedges and HW
        consolidateDirectedHG();
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
//
//		System.out.printf("sum of Htail[198,:]: %.6f\n", Htail.getRow(198).fold(Vectors.asSumAccumulator(0)));
//		for (int e = 0; e < E; ++e) {
//			double v = Htail.get(198, e);
//			if (v != 0)
//				System.out.printf("Htail[198,%d] = %.16f\n", e,v);
//		}
//		System.out.printf("OutDv(198) = %.6f\n", getTailVertexDegree(198,false));
//		List<Integer> heids = getIncidentOutHyperedgeIDs(198);
//		for (int heid : heids) {
//			System.out.printf("HWtail[198,%d] = %.16f\n", heid, HWtail.get(198, heid));
//		}
    }


    /**
     * Given H, WE, WL, hyperedgeWeighted and vertexWeighted either explicitly or in default
     * calculate hyperedges and HW
     */
    private void consolidateDirectedHG() {	// 需要修改，考虑weighting
        if (hyperedges != null)
            hyperedges.clear();
        hyperedges = new ArrayList<DirectedHyperEdge>(E);
        for (int e = 0; e < E; ++e) {
            DirectedHyperEdge he = new DirectedHyperEdge(e);
            hyperedges.add(he);
        }
        //
        MatrixIterator itTail = null;
        if (HWtail instanceof DenseMatrix) {
            itTail = HWtail.iterator();
        } else {
            if (HWtail instanceof RowMajorSparseMatrix)
                itTail = ((CRSMatrix) HWtail).nonZeroRowMajorIterator();
            else
                itTail = ((CCSMatrix) HWtail).nonZeroColumnMajorIterator();
        }
        while (itTail.hasNext()) {
            double w = itTail.next();
            int v = itTail.rowIndex();
            int e = itTail.columnIndex();
//			if (v == 198 && e == 62540)
//				System.out.println("v = " + v + "; e = " + e);
            if (w > 0)
                hyperedges.get(e).addTailVertex(v);
        }
        //
        MatrixIterator itHead = null;
        if (HWhead instanceof DenseMatrix) {
            itHead = HWhead.iterator();
        } else {
            if (HWhead instanceof RowMajorSparseMatrix)
                itHead = ((CRSMatrix) HWhead).nonZeroRowMajorIterator();
            else
                itHead = ((CCSMatrix) HWhead).nonZeroColumnMajorIterator();
        }
        while (itHead.hasNext()) {
            double w = itHead.next();
            int v = itHead.rowIndex();
            int e = itHead.columnIndex();
            if (w > 0)
                hyperedges.get(e).addHeadVertex(v);
        }
        //
//        for (int e = 0; e < E; ++e) {
//            hyperedges.get(e).consolidate();
//        }
    }


    @Override
    public final Matrix getOutAdjacencyMatrix(boolean vertexWeighted) {
        if (vertexWeighted)
            return HWtail;
        else
            return Htail;
    }
    @Override
    public final Matrix getInAdjacencyMatrix(boolean vertexWeighted) {
        if (vertexWeighted)
            return HWhead;
        else
            return Hhead;
    }
    /*
        public boolean isHyperedgeWeighted() {
            return hyperedgeWeighted;
        }
    */
    public final Vector getHyperedgeWeights() {
        return WE;
    }

    public double getHyperedgeWeight(int e) {
        return WE.get(e);
    }

    public double getHyperedgeWeight(int e, boolean weighted) {
        if (weighted)
            return WE.get(e);
        else
            return 1;
    }
/*
	public boolean isVertexWeighted() {
		return vertexWeighted;
	}
*/
    /**
     * Equivalent to getOutAdjacencyMatrix().
     * @return
     */
    public final Matrix getTailVertexWeights() {
        return HWtail;
    }

    /**
     * Equivalent to getInAdjacencyMatrix().
     * @return
     */
    public final Matrix getHeadVertexWeights() {
        return HWhead;
    }

    /**
     * Get the weights of the v-th tail vertex corresponding to all the hyperedges
     * @param v
     * @return a vector of weights of the v-th tail vertex corresponding to all the hyperedges
     */
    public final Vector getTailVertexWeights(int v) {
        return HWtail.getRow(v);
    }

    /**
     * Get the weights of the v-th head vertex corresponding to all the hyperedges
     * @param v
     * @return a vector of weights of the v-th head vertex corresponding to all the hyperedges
     */
    public final Vector getHeadVertexWeights(int v) {
        return HWhead.getRow(v);
    }

    /**
     * Get the weight of the v-th tail vertex corresponding to the e-th hyperedge
     * @param v
     * @param e
     * @return the weight of the v-th tail vertex corresponding to the e-th hyperedge
     */
    public double getTailVertexWeight(int v, int e) {
        return HWtail.get(v, e);
    }

    public double getTailVertexWeight(int v, int e, boolean vertexWeighted) {
        if (vertexWeighted)
            return HWtail.get(v, e);
        else
            return Htail.get(v, e);
    }

    /**
     * Get the weight of the v-th head vertex corresponding to the e-th hyperedge
     * @param v
     * @param e
     * @return the weight of the v-th head corresponding to the e-th hyperedge
     */
    public double getHeadVertexWeight(int v, int e) {
        return HWhead.get(v, e);
    }

    public double getHeadVertexWeight(int v, int e, boolean vertexWeighted) {
        if (vertexWeighted)
            return HWhead.get(v, e);
        else
            return Hhead.get(v, e);
    }


    @Override
    public List<Double> getHyperedgeOutDegrees(boolean vertexWeighted) {
        List<Double> vexEdgeDeg = new ArrayList<Double>(E);
        for (int e = 0; e < E; ++e) {
            vexEdgeDeg.add( getHyperedgeOutDegree(e, vertexWeighted) );
        }
        return vexEdgeDeg;
    }

    /**
     * Get the degree of the e-th hyperedge, according to the parameter weighted
     * deg(e) = SIGMA_{v incident to e}{w(v,e)}    w(v,e) = HW(v,e) or H(v,e)
     * @param e
     * @param vertexWeighted if weighted, w(v,e) = HW(v,e); otherwise w(v,e) = HW(v,e)
     * @return
     */
    @Override
    public double getHyperedgeOutDegree(int e, boolean vertexWeighted) {
        double deg = 0;
        for (int v : hyperedges.get(e).getTailVertexIDs()) {
            if (vertexWeighted)
                deg += HWtail.get(v, e);
            else
                deg++;
        }
        return deg;
    }

    @Override
    public List<Double> getHyperedgeInDegrees(boolean vertexWeighted) {
        List<Double> vexEdgeDeg = new ArrayList<Double>(E);
        for (int e = 0; e < E; ++e) {
            vexEdgeDeg.add( getHyperedgeInDegree(e, vertexWeighted) );
        }
        return vexEdgeDeg;
    }

    /**
     * Get the degree of the e-th hyperedge, according to the parameter weighted
     * deg(e) = SIGMA_{v incident to e}{w(v,e)}    w(v,e) = HW(v,e) or H(v,e)
     * @param e
     * @param vertexWeighted if weighted, w(v,e) = HW(v,e); otherwise w(v,e) = HW(v,e)
     * @return
     */
    @Override
    public double getHyperedgeInDegree(int e, boolean vertexWeighted) {
        double deg = 0;
        for (int v : hyperedges.get(e).getHeadVertexIDs()) {
            if (vertexWeighted)
                deg += HWhead.get(v, e);
            else
                deg++;
        }
        return deg;
    }

    @Override
    public List<Double> getHyperedgeDegrees(boolean vertexWeighted) {
        List<Double> degs = new ArrayList<Double>(E);
        for (int e = 0; e < E; ++e) {
            double deg = getHyperedgeDegree(e, vertexWeighted);
            degs.add(deg);
        }
        return degs;
    }

    @Override
    public double getHyperedgeDegree(int e, boolean vertexWeighted) {
        double inDeg = getHyperedgeInDegree(e, vertexWeighted);
        double outDeg = getHyperedgeOutDegree(e, vertexWeighted);
        double deg = inDeg + outDeg;
        return deg;
    }


    @Override
    public List<Double> getTailVertexDegrees(boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        List<Double> degs = new ArrayList<Double>(V);
        for (int v = 0; v < V; ++v) {
            double deg = getTailVertexDegree(v, hyperedgeWeighted);
            degs.add(deg);
        }
        return degs;
    }
    /**
     *
     * @param v the id of vertex which is as a tail vertex
     * @return
     */
    @Override
    public double getTailVertexDegree(int v, boolean hyperedgeWeighted) {
        double deg = 0;
        List<Integer> heids = getIncidentOutHyperedgeIDs(v);
        for (int heid : heids) {
            if (hyperedgeWeighted)
                deg += WE.get(heid);
            else
                deg++; //+= H.get(v, heid);
        }
//		deg = heids.size();
//		if (v == 198)
//			System.out.println(deg);
        return deg;
    }


    @Override
    public List<Double> getHeadVertexDegrees(boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        List<Double> degs = new ArrayList<Double>(V);
        for (int v = 0; v < V; ++v) {
            double deg = getHeadVertexDegree(v, hyperedgeWeighted);
            degs.add(deg);
        }
        return degs;
    }
    /**
     *
     * @param v the id of vertex which is as a head vertex
     * @return
     */
    @Override
    public double getHeadVertexDegree(int v, boolean hyperedgeWeighted) {
        double deg = 0;
        List<Integer> heids = getIncidentInHyperedgeIDs(v);
        for (int heid : heids) {
            if (hyperedgeWeighted)
                deg += WE.get(heid);
            else
                deg++; //+= H.get(v, heid);
        }
        return deg;
    }


    @Override
    public List<Double> getVertexDegrees(boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        List<Double> degs = new ArrayList<Double>(V);
        for (int v = 0; v < V; ++v) {
            double deg = getVertexDegree(v, hyperedgeWeighted);
            degs.add(deg);
        }
        return degs;
    }

    @Override
    public double getVertexDegree(int v, boolean hyperedgeWeighted) {
        // TODO Auto-generated method stub
        double tailDeg = getTailVertexDegree(v, hyperedgeWeighted);
        double headDeg = getHeadVertexDegree(v, hyperedgeWeighted);
        double deg = tailDeg + headDeg;
        return deg;
    }


    /**
     * Get the ids of the hyperedges for which vertex v is a tail vertex, i.e. these hyperedges are out edges from v
     * @param v
     * @return
     */
    public final List<Integer> getIncidentOutHyperedgeIDs(int v) {
        List<Integer> heids = new ArrayList<Integer>();
        VectorIterator it;
        Vector row = Htail.getRow(v);
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
     * Get the ids of the hyperedges for which vertex v is a tail vertex
     * @param v
     * @return
     */
    public final List<DirectedHyperEdge> getIncidentOutHyperedges(int v) {
        List<Integer> heids = getIncidentOutHyperedgeIDs(v);
        List<DirectedHyperEdge> hes = new ArrayList<DirectedHyperEdge>(heids.size());
        for (int heid : heids) {
            hes.add(hyperedges.get(heid));
        }
        return hes;
    }

    /**
     * Get the ids of the hyperedges for which vertex v is a head vertex, i.e. these hyperedges are in edges to v
     * @param v
     * @return
     */
    public final List<Integer> getIncidentInHyperedgeIDs(int v) {
        List<Integer> heids = new ArrayList<Integer>();
        VectorIterator it = Hhead.iteratorOfRow(v);
        Vector row = Hhead.getRow(v);
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
     * Get the ids of the hyperedges for which vertex v is a head vertex
     * @param v
     * @return
     */
    public final List<DirectedHyperEdge> getIncidentInHyperedges(int v) {
        List<Integer> heids = getIncidentInHyperedgeIDs(v);
        List<DirectedHyperEdge> hes = new ArrayList<DirectedHyperEdge>(heids.size());
        for (int heid : heids) {
            hes.add(hyperedges.get(heid));
        }
        return hes;
    }

    /**
     * Get the ids of the tail vertices incident to the e-th hyperedge
     * @param e
     * @return
     */
    public final List<Integer> getIncidentTailVertices(int e) {
        return hyperedges.get(e).getTailVertexIDs();
    }

    /**
     * Get the ids of the head vertices incident to the e-th hyperedge
     * @param e
     * @return
     */
    public final List<Integer> getIncidentHeadVertices(int e) {
        return hyperedges.get(e).getHeadVertexIDs();
    }




    /**
     * Find sink tail or head vertices according to input parameter <i>tailOrHead</i>
     * Sink tail vertex: not incident to any hyperedge as tail vertex
     * Sink head vertex: not incident to any hyperedge as head vertex
     * @param tailOrHead true if finding sink tail vertices, false if finding sink head vertices
     * @return
     */
    public Vector findSinkVertices(boolean tailOrHead) {
        Vector sinks = BasicVector.zero(V);
        for (int v = 0; v < V; ++v) {
            if (isSinkVertex(v, tailOrHead))
                sinks.set(v, 1);
        }
        return sinks;
    }

    public boolean isSinkVertex(int v, boolean tailOrHead) {
        if (tailOrHead) {
            double dv = getTailVertexDegree(v, hyperedgeWeighted);	// number of directed hyperedges out of vertex v
            return Utility.isZero(dv);
        } else {
            double dv = getHeadVertexDegree(v, hyperedgeWeighted);	// number of directed hyperedges out of vertex v
            return Utility.isZero(dv);
        }
    }

    public Vector findSinkHyperedges(boolean inOrOut) {
        Vector sinks = BasicVector.zero(E);
        for (int e = 0; e < E; ++e) {
            if (isSinkHyperedge(e, inOrOut))
                sinks.set(e, 1);
        }
        return sinks;
    }

    public boolean isSinkHyperedge(int e, boolean inOrOut) {
        if (inOrOut) {
            double de = getHyperedgeInDegree(e, vertexWeighted);	// number of tail vertices incident to hyperedge e
            return Utility.isZero(de);
        } else {
            double de = getHyperedgeOutDegree(e, vertexWeighted);	// number of head vertices incident to hyperedge e
            return Utility.isZero(de);
        }
    }



    /**
     * Normalize adjacency matrix HW to transition matrix P, actually also H2A
     *
     * 1. neither hyperedge-weighted nor vertex-weighted <br>
     * P = Dvout^(-1) * Htail * Dein^(-1) * Hhead^(T) <br>
     * where <br>
     * &nbsp;&nbsp;Htail[u,e] = 1 if e starts from u; Hhead[u,e] = 1 if e ends with u <br>
     * &nbsp;&nbsp;Dvout[u] = SIGMA_e{ Htail[u,e] } = |E(u)| <br>
     * &nbsp;&nbsp;Dein[e] = SIGMA_v{ Hhead[v,e] } = |V(e)| <br>
     * &nbsp;&nbsp;We = diag(WE) = I <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;WE = constant 1-vector <br>
     * &nbsp;&nbsp;HW = HE <br>
     *
     * 2. hyperedge-weighted but not vertex weighted <br>
     * P = Dvout^(-1) * Htail * We * Dein^(-1) * Hhead^(T) <br>
     * where <br>
     * &nbsp;&nbsp;Htail[u,e] = 1 if e starts from u; Hhead[u,e] = 1 if e ends with u <br>
     * &nbsp;&nbsp;Dvout[u] = SIGMA_e{ WE[e] * Htail[u,e] } <br>
     * &nbsp;&nbsp;Dein[e] = SIGMA_v{ Hhead[v,e] } = |V(e)| <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;HW = HE <br>
     *
     * 3. both hyperedge-weighted and vertex-weighted <br>
     * P = Dvout^(-1) * Htail * We * Dein^(-1) * HWhead^(T) <br>
     * where <br>
     * &nbsp;&nbsp;Htail[u,e] = 1 if e starts from u; Hhead[u,e] = 1 if e ends with u <br>
     * &nbsp;&nbsp;Dvout[u] = SIGMA_e{ WE[e] * Htail[u,e] } <br>
     * &nbsp;&nbsp;Dein[e] = SIGMA_v{ HWhead[v,e] }  <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;and <br>
     * &nbsp;&nbsp;HWtail[u,e] > 0, HWhead[u,e] > 0, Htail[u,e] = 1 if e starts from u; Hhead[u,e] = 1 if e ends with u <br>
     *
     * 4. All in one, the following is a general form <br>
     * P = Dvout^(-1) * Htail * We * Dein^(-1) * HWhead^(T) <br>
     * where <br>
     * &nbsp;&nbsp;Htail[u,e] = 1 if e starts from u; Hhead[u,e] = 1 if e ends with u <br>
     * &nbsp;&nbsp;Dvout[u] = SIGMA_e{ WE[e] * Htail[u,e] } <br>
     * &nbsp;&nbsp;Dein[e] = SIGMA_v{ HWhead[v,e] } <br>
     * &nbsp;&nbsp;We = diag(WE) <br>
     * &nbsp;&nbsp;HWtail[u,e] > 0, HWhead[u,e] > 0, Htail[u,e] = 1 if e starts from u; Hhead[u,e] = 1 if e ends with u <br>
     * &nbsp;If the hypergraph is neither hyperedge-weighted nor vertex-weighted <br>
     * &nbsp;&nbsp;Dvout[u] = |E(u)| <br>
     * &nbsp;&nbsp;Dein[e] = |V(e)| <br>
     * &nbsp;&nbsp;We = I <br>
     * &nbsp;&nbsp;HW = H <br>
     * &nbsp;If the hypergraph is hyperedge-weighted but not vertex-weighted <br>
     * &nbsp;&nbsp;Dein[e] = |V(e)| <br>
     * &nbsp;&nbsp;HWtail = Htail, HWhead = Hhead <br>
     * &nbsp;If the hypergraph is both hyperedge-weighted and vertex-weighted <br>
     * @return
     */
    @Override
    public Matrix normalize(boolean hyperedgeWeighted, boolean vertexWeighted) {
        if (allowSelfLoop)
            return normalizeWithSelfLoop(hyperedgeWeighted, vertexWeighted);
        else
            return normalizeWithoutSelfLoop(hyperedgeWeighted, vertexWeighted);
    }

    private Matrix normalizeWithSelfLoop(boolean hyperedgeWeighted, boolean vertexWeighted) {
        Matrix P = null;
//		int V = this.getNumberOfVertices(), E = this.getNumberOfHyperedges();
        // Htail, Hhead and HWtail, HWhead
        Matrix Htail = this.getOutAdjacencyMatrix(false);
        Matrix HWhead = this.getInAdjacencyMatrix(vertexWeighted);
        boolean vertexIsTailOrHead = true;
        Vector sinkVertices = this.findSinkVertices(vertexIsTailOrHead);		// sink tail vertices
        boolean hyperedgeisInOrOut = true;
        Vector sinkHyperedges = this.findSinkHyperedges(hyperedgeisInOrOut);		// sink in hyperedges
//		// Inverse of Dvout
//		Matrix DvoutI = CRSMatrix.diagonal(V, 1);
//		for (int v = 0; v < V; ++v) {
//			if (Utility.isZero(sinks.get(v))) {	// not sink
//				double dv = this.getTailVertexDegree(v, hyperedgeWeighted);
//				DvoutI.set(v, v, 1.0/dv);
//			} else
//				DvoutI.set(v, v, 0);
//		}
//		// We: diagonal matrix of hyperedge weights
//		Matrix We = CRSMatrix.diagonal(E, 1);
//		for (int e = 0; e < E; ++e) {
//			double we = this.getHyperedgeWeight(e);
//			We.set(e, e, we);
//		}
//		// Inverse of Dein
//		Matrix DeinI = CRSMatrix.diagonal(E, 1);
//		for (int e = 0; e < E; ++e) {
//			double de = this.getHyperedgeInDegree(e, hyperedgeWeighted);
//			DeinI.set(e, e, 1.0/de);
//		}
//		//
//		P = (DvoutI.multiply(Htail).multiply(We)).multiply(DeinI.multiply(HWhead.transpose()));
        // reduce the number of multiplications: codes ...
        Matrix P1 = Htail.copy();
        Vector OutDv = BasicVector.zero(V);
        for (int v = 0; v < V; ++v) {
            if (Utility.isZero(sinkVertices.get(v))) {	// not sink vertices
                double dv = this.getTailVertexDegree(v, hyperedgeWeighted);
                OutDv.set(v, dv);
            }
        }
        Vector We = BasicVector.constant(E, 1);
        Vector InDe = BasicVector.zero(E);
        for (int e = 0; e < E; ++e) {
            if (hyperedgeWeighted) {
                double we = this.getHyperedgeWeight(e);
                We.set(e, we);
            }
            if (Utility.isZero(sinkHyperedges.get(e))) {	// not sink hyperedges
                double de = this.getHyperedgeInDegree(e, vertexWeighted);
                InDe.set(e, de);
            }
        }
        // P = OutDv^(-1) * Htail * We * InDe^(-1) * HWhead^(T)
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
        // 1. OutDv^(-1) * Htail * We * InDe^(-1);  P1 = Htail
        while (itP1.hasNext()) {
            double wv = itP1.next();
            int v = itP1.rowIndex();
            if (!Utility.isZero(sinkVertices.get(v)))	// sink vertices
                continue;
            wv = wv / OutDv.get(v);				// divide each v-th row of P1 by OutDv(v)
            int e = itP1.columnIndex();
            if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
                continue;
            wv = wv * We.get(e) / InDe.get(e);			// multiply each e-th cloumn by We(e)
            itP1.set(wv);
        }
//		for (int i = 0; i < 500; ++i) {
//			System.out.printf("sum of P1[%d,:]: %f\n", i, P1.getRow(i).fold(Vectors.asSumAccumulator(0)));
//		}
//		System.out.println();
        //
//		Matrix P2 = HWhead.copy();
//		MatrixIterator itP2 = null;
//		if (P2 instanceof DenseMatrix) {
//			itP2 = P2.iterator();
//		} else {
//			if (P2 instanceof RowMajorSparseMatrix) {
//				itP2 = ((CRSMatrix) P2).nonZeroRowMajorIterator();
//			} else {
//				itP2 = ((CCSMatrix) P2).nonZeroColumnMajorIterator();
//			}
//		}
//		while (itP2.hasNext()) {
//			double wv = itP2.next();
//			int e = itP2.columnIndex();
//			if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
//				continue;
//			wv = wv / InDe.get(e);	// multiply each e-th cloumn by We(e) and divide the same column by InDe(e)
//			itP2.set(wv);
//		}
//		for (int i = 0; i < 500; ++i) {
//			System.out.printf("sum of P2[:,%d]: %f\n", i, P2.getColumn(i).fold(Vectors.asSumAccumulator(0)));
//		}
//		System.out.println();
        // 2. OutDv^(-1) * Htail * We * InDe^(-1) * HWhead^(T);
        //
        P = P1.multiply(HWhead.transpose());
//		for (int i = 0; i < 200; ++i) {
//			System.out.printf("sum of row[%d]: %f\n", i, P.getRow(i).fold(Vectors.asSumAccumulator(0)));
//		}
//		System.out.println();
        return P;
    }

    private Matrix normalizeWithoutSelfLoop(boolean hyperedgeWeighted, boolean vertexWeighted) {
        Matrix P = null;
        //src_vid --> list of pairs of <trg_vid, trans_prob>
        //             trg_vid  trans_prob
        Map<Pair<Integer, Integer>, RealAccumulator> TransProbsMap = new HashMap<>();
        for (int src = 0; src < V; ++src) { // transition probabilities from source node src to all the "incident" target nodes
            // From tail vertex src to its incident out hyperedges
            List<Integer> outHyperedgeIDs = this.getIncidentOutHyperedgeIDs(src);
            double sumOutHyperedgeWeight = this.getTailVertexDegree(src, hyperedgeWeighted);
            int outE = outHyperedgeIDs.size();
            for (int outHeid : outHyperedgeIDs) { // transition probabilities from source node src via each out hypergraph
                double outHyperedgeWeight = this.getHyperedgeWeight(outHeid, hyperedgeWeighted);
                // calculate transition probability from source node src to this out hyperedge outHeid
                double transProbSrc2HE = outHyperedgeWeight / sumOutHyperedgeWeight;
                // From hyperedge outHeid to its incident head vertices
                List<Integer> headVertexIDs = this.getIncidentHeadVertices(outHeid);
                headVertexIDs.remove(new Integer(src));		// self-loop not allowed
                if (headVertexIDs.isEmpty()) continue;		// no target nodes after removeing self-loop
                double sumHeadVertexWeight = 0;	// calculate indegree of this hyperedge without self-loop
                for (int trg : headVertexIDs) {
                    double headVertexWeight = this.getHeadVertexWeight(trg, outHeid, vertexWeighted);
                    sumHeadVertexWeight += headVertexWeight;
                }
                int outV = headVertexIDs.size();
                for (int trg : headVertexIDs) { // transition probabilities to target nodes
                    // calculate transition probability from hyperedge outHeid to target node trg
                    double headVertexWeight = this.getHeadVertexWeight(trg, outHeid, vertexWeighted);
                    double transProbHE2Trg = headVertexWeight / sumHeadVertexWeight;
                    // calculate transition probability from source node src to target node trg
                    double transProbSrc2Trg = transProbSrc2HE * transProbHE2Trg;
                    // Accumulate this portion of transition probability in TransProbsMap
                    Pair<Integer, Integer> st = new Pair<Integer, Integer>(src, trg);
                    if (!TransProbsMap.containsKey(st)) {
                        TransProbsMap.put(st, new RealAccumulator(transProbSrc2Trg));
                    } else {
                        TransProbsMap.get(st).add(transProbSrc2Trg);
                    }
                } // for each trg
            } // for each outHeid
        } // for each src
        // TransProbsMap --> TransProbsList
        int cardP = TransProbsMap.size();
//		List<Pair<Pair<Integer, Integer>, Double>> TransProbsList = new ArrayList<Pair<Pair<Integer, Integer>, Double>>(cardP);
//		for (Pair<Integer, Integer> st : TransProbsMap.keySet()) {
//			Pair<Pair<Integer, Integer>, Double> triple =
//					new Pair<Pair<Integer, Integer>, Double>(st, TransProbsMap.get(st).getRealValue());
//			TransProbsList.add(triple);
//		}
//		Collections.sort(TransProbsList);
        List<List<Pair<Integer, Double>>> TransProbsList = new ArrayList<List<Pair<Integer, Double>>>(V);
        for (int src = 0; src < V; ++src) {
            TransProbsList.add( new LinkedList<Pair<Integer, Double>>() );
        }
        for (Pair<Integer, Integer> st : TransProbsMap.keySet()) {
            Pair<Pair<Integer, Integer>, Double> triple =
                    new Pair<Pair<Integer, Integer>, Double>(st, TransProbsMap.get(st).getRealValue());
            int src = triple.getValue0().getValue0();
            int trg = triple.getValue0().getValue1();
            double prob = triple.getValue1();
            TransProbsList.get(src).add( new Pair<Integer, Double>(trg, prob) );
        }
        for (int src = 0; src < V; ++src) {
            Collections.sort(TransProbsList.get(src));
        }
        // TransProbsList --> P
        int[] columnIndices = new int[cardP];
        double[] values = new double[cardP];
        int[] rowPointers = new int[V+1];
        rowPointers[0] = 0;
//		int preSrc = 0;
//		for (int i = 0; i < cardP; ++i) {
//			Pair<Pair<Integer, Integer>, Double> triple = TransProbsList.get(i);
//			int src = triple.getValue0().getValue0();
//			int trg = triple.getValue0().getValue1();
//			double prob = triple.getValue1();
//			columnIndices[i] = trg;
//			values[i] = prob;
//			if (src != preSrc)
//		}
        int k = 0;
        for (int src = 0; src < V; ++src) {
            List<Pair<Integer, Double>> outLinks = TransProbsList.get(src);
            rowPointers[src+1] = rowPointers[src] + outLinks.size();
            for (Pair<Integer, Double> outLink : outLinks) {
                int trg = outLink.getValue0();
                double prob = outLink.getValue1();
                columnIndices[k] = trg;
                values[k++] = prob;
            }
        }
        P = new CRSMatrix(V, V, cardP, values, columnIndices, rowPointers);
        //
        return P;
    }


    // need codes ...  20150904
    /**
     *
     * @param H2A  whether normalize it as a hub-to-authority transition matrix as PageRank does
     * or as an authori-to-hub transition matrix as HITS does
     * @return
     */
    @Override
    public Matrix normalize(boolean H2A, boolean hyperedgeWeighted, boolean vertexWeighted) {
        // TODO Auto-generated method stub
        if (H2A)	// H2A : hub-to-authority transition matrix
            return normalize(hyperedgeWeighted, vertexWeighted);
        else
            return normalizeA2H(hyperedgeWeighted, vertexWeighted);
    }



    private Matrix normalizeA2H(boolean hyperedgeWeighted, boolean vertexWeighted) {
        if (allowSelfLoop)
            return normalizeWithSelfLoopA2H(hyperedgeWeighted, vertexWeighted);
        else
            return normalizeWithoutSelfLoopA2H(hyperedgeWeighted, vertexWeighted);

    }

    private Matrix normalizeWithSelfLoopA2H(boolean hyperedgeWeighted, boolean vertexWeighted) {
        // A2H : authority-to-hub transition matrix
        // NormH2A = OutDv^(-1) * Htail * We * InDe^(-1) * HWhead^(T);
        // while NormA2H = InDv^(-1) * Hhead * We * OutDe^(-1) * HWtail^(T);
        Matrix P = null;
        Matrix Hhead = this.getInAdjacencyMatrix(false);
        Matrix HWtail = this.getOutAdjacencyMatrix(vertexWeighted);
        boolean vertexisTailOrHead = false;
        Vector sinkVertices = this.findSinkVertices(vertexisTailOrHead);		// sink tail vertices
        boolean hyperedgeisInOrOut = false;
        Vector sinkHyperedges = this.findSinkHyperedges(hyperedgeisInOrOut);		// sink in hyperedges
        // reduce the number of multiplications: codes ...
        // NormA2H = InDv^(-1) * Hhead * We * OutDe^(-1) * HWtail^(T);
        Matrix P1 = Hhead.copy();
        Vector InDv = BasicVector.zero(V);
        for (int v = 0; v < V; ++v) {
            if (Utility.isZero(sinkVertices.get(v))) {	// not sink vertices
                double dv = this.getHeadVertexDegree(v, hyperedgeWeighted);
                InDv.set(v, dv);
            }
        }
        Vector We = BasicVector.constant(E, 1);
        Vector OutDe = BasicVector.zero(E);
        for (int e = 0; e < E; ++e) {
            if (hyperedgeWeighted) {
                double we = this.getHyperedgeWeight(e);
                We.set(e, we);
            }
            if (Utility.isZero(sinkHyperedges.get(e))) {	// not sink hyperedges
                double de = this.getHyperedgeOutDegree(e, vertexWeighted);
                OutDe.set(e, de);
            }
        }
        // P = InDv^(-1) * Hhead * We * OutDe^(-1) * HWtail^(T)
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
        // 1. InDv^(-1) * Hhead * We * OutDe^(-1);  P1 = Hhead
        while (itP1.hasNext()) {
            double wv = itP1.next();
            int v = itP1.rowIndex();
            if (!Utility.isZero(sinkVertices.get(v)))	// sink vertices
                continue;
            wv = wv / InDv.get(v);				// divide each v-th row of P1 by OutDv(v)
            int e = itP1.columnIndex();
            if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
                continue;
            wv = wv * We.get(e) / OutDe.get(e);	// multiply each e-th cloumn by We(e) and divide the same column by InDe(e)
            itP1.set(wv);
        }
//		for (int i = 0; i < 200; ++i) {
//			System.out.printf("sum of P1[%d,:]: %f\n", i, P1.getRow(i).fold(Vectors.asSumAccumulator(0)));
//		}
        //
//		Matrix P2 = HWtail.copy();
//		MatrixIterator itP2 = null;
//		if (P2 instanceof DenseMatrix) {
//			itP2 = P2.iterator();
//		} else {
//			if (P2 instanceof RowMajorSparseMatrix) {
//				itP2 = ((CRSMatrix) P2).nonZeroRowMajorIterator();
//			} else {
//				itP2 = ((CCSMatrix) P2).nonZeroColumnMajorIterator();
//			}
//		}
//		while (itP2.hasNext()) {
//			double wv = itP2.next();
//			int e = itP2.columnIndex();
//			if (!Utility.isZero(sinkHyperedges.get(e)))	// sink hyperedges
//				continue;
//			wv = wv / OutDe.get(e);	// multiply each e-th cloumn by We(e) and divide the same column by InDe(e)
//			itP2.set(wv);
//		}
//		for (int i = 0; i < 200; ++i) {
//			System.out.printf("sum of P2[:,%d]: %f\n", i, P2.getColumn(i).fold(Vectors.asSumAccumulator(0)));
//		}
        // 2. InDv^(-1) * Hhead * We * OutDe^(-1) * HWtail^(T)
        P = P1.multiply(HWtail.transpose());
//		P = P1.multiply(P2.transpose());
//		MatrixUtil.pprint(P);
//		for (int i = 0; i < 200; ++i) {
//			System.out.printf("sum of row[%d]: %f\n", i, P.getRow(i).fold(Vectors.asSumAccumulator(0)));
//		}
        return P;
    }

    // reverse citation hypernetwork
    private Matrix normalizeWithoutSelfLoopA2H(boolean hyperedgeWeighted, boolean vertexWeighted) {
        Matrix P = null;
        //src_vid --> list of pairs of <trg_vid, trans_prob>
        //             trg_vid  trans_prob
        Map<Pair<Integer, Integer>, RealAccumulator> TransProbsMap = new HashMap<Pair<Integer, Integer>, RealAccumulator>();
        for (int src = 0; src < V; ++src) { // transition probabilities from source node src to all the "incident" target nodes
            // From head vertex src to its incident in hyperedges
            List<Integer> inHyperedgeIDs = this.getIncidentInHyperedgeIDs(src);
            double sumInHyperedgeWeight = this.getHeadVertexDegree(src, hyperedgeWeighted);
            int inE = inHyperedgeIDs.size();
            for (int inHeid : inHyperedgeIDs) { // transition probabilities from source node src via each in hyperedge
                double inHyperedgeWeight = this.getHyperedgeWeight(inHeid, hyperedgeWeighted);
                // calculate transition probability from source node src to this in hyperedge inHeid
                double transProbSrc2HE = inHyperedgeWeight / sumInHyperedgeWeight;
                // From hyperedge inHeid to its incident tail vertices
                List<Integer> tailVertexIDs = this.getIncidentTailVertices(inHeid);

                tailVertexIDs.remove(new Integer(src));		// self-loop not allowed
                if (tailVertexIDs.isEmpty()) continue;		// no target nodes after removeing self-loop
                double sumTailVertexWeight = 0;	// calculate indegree of this hyperedge without self-loop
                for (int trg : tailVertexIDs) {
                    double tailVertexWeight = this.getTailVertexWeight(trg, inHeid, vertexWeighted);
                    sumTailVertexWeight += tailVertexWeight;
                }
                int inV = tailVertexIDs.size();
                for (int trg : tailVertexIDs) { // transition probabilities to target nodes
                    // calculate transition probability from hyperedge inHeid to target node trg
                    double tailVertexWeight = this.getTailVertexWeight(trg, inHeid, vertexWeighted);
                    double transProbHE2Trg = tailVertexWeight / sumTailVertexWeight;
                    // calculate transition probability from source node src to target node trg
                    double transProbSrc2Trg = transProbSrc2HE * transProbHE2Trg;
                    // Accumulate this portion of transition probability in TransProbsMap
                    Pair<Integer, Integer> st = new Pair<Integer, Integer>(src, trg);
                    if (!TransProbsMap.containsKey(st)) {
                        TransProbsMap.put(st, new RealAccumulator(transProbSrc2Trg));
                    } else {
                        TransProbsMap.get(st).add(transProbSrc2Trg);
                    }
                } // for each trg
            } // for each outHeid
        } // for each src
        // TransProbsMap --> TransProbsList
        int cardP = TransProbsMap.size();
//		List<Pair<Pair<Integer, Integer>, Double>> TransProbsList = new ArrayList<Pair<Pair<Integer, Integer>, Double>>(cardP);
//		for (Pair<Integer, Integer> st : TransProbsMap.keySet()) {
//			Pair<Pair<Integer, Integer>, Double> triple =
//					new Pair<Pair<Integer, Integer>, Double>(st, TransProbsMap.get(st).getRealValue());
//			TransProbsList.add(triple);
//		}
//		Collections.sort(TransProbsList);
        List<List<Pair<Integer, Double>>> TransProbsList = new ArrayList<List<Pair<Integer, Double>>>(V);
        for (int src = 0; src < V; ++src) {
            TransProbsList.add( new LinkedList<Pair<Integer, Double>>() );
        }
        for (Pair<Integer, Integer> st : TransProbsMap.keySet()) {
            Pair<Pair<Integer, Integer>, Double> triple =
                    new Pair<Pair<Integer, Integer>, Double>(st, TransProbsMap.get(st).getRealValue());
            int src = triple.getValue0().getValue0();
            int trg = triple.getValue0().getValue1();
            double prob = triple.getValue1();
            TransProbsList.get(src).add( new Pair<Integer, Double>(trg, prob) );
        }
        for (int src = 0; src < V; ++src) {
            Collections.sort(TransProbsList.get(src));
        }
        // TransProbsList --> P
        int[] columnIndices = new int[cardP];
        double[] values = new double[cardP];
        int[] rowPointers = new int[V+1];
        rowPointers[0] = 0;
//		int preSrc = 0;
//		for (int i = 0; i < cardP; ++i) {
//			Pair<Pair<Integer, Integer>, Double> triple = TransProbsList.get(i);
//			int src = triple.getValue0().getValue0();
//			int trg = triple.getValue0().getValue1();
//			double prob = triple.getValue1();
//			columnIndices[i] = trg;
//			values[i] = prob;
//			if (src != preSrc)
//		}
        int k = 0;
        for (int src = 0; src < V; ++src) {
            List<Pair<Integer, Double>> outLinks = TransProbsList.get(src);
            rowPointers[src+1] = rowPointers[src] + outLinks.size();
            for (Pair<Integer, Double> outLink : outLinks) {
                int trg = outLink.getValue0();
                double prob = outLink.getValue1();
                columnIndices[k] = trg;
                values[k++] = prob;
            }
        }
        P = new CRSMatrix(V, V, cardP, values, columnIndices, rowPointers);
        //
        return P;
    }
}
