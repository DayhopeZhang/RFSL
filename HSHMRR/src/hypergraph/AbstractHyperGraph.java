// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.List;

import org.la4j.Matrix;
import org.la4j.Vector;

public abstract class AbstractHyperGraph {

    protected int V, E;		// number of Vertices and Edges
    protected Vector WE;	// weights of hyperedges
    // IF hyperedgeWeighted assigned to FALSE, then set WE to a unit vector (1-vector)
    // OTHERWISE
    //   IF constructor assigns WE, the use this WE as hyperedge weight
    //   OTHERWISE use HW to calculate hyperedge weight, that is hyperedgeWeight
    protected boolean hyperedgeWeighted;
    protected boolean vertexWeighted;

    public boolean isVertexWeighted() {
        return vertexWeighted;
    }

    public boolean isHyperedgeWeighted() {
        return hyperedgeWeighted;
    }

    public int getNumberOfHyperedges() {
        return E;
    }

    public int getNumberOfVertices() {
        return V;
    }


    ///////////////////////////////////////////////////////////////////////////
    //////// methods for calculating hyperedge degrees
    //////////////////////////////////////////////////////////////////////////
    // hyperedge in degrees
    public final List<Double> getHyperedgeInDegrees() {
        return getHyperedgeInDegrees(vertexWeighted);
    }
    abstract public List<Double> getHyperedgeInDegrees(boolean vertexWeighted);
    // hyperedge in degree
    public double getHyperedgeInDegree(int e) {
        return getHyperedgeInDegree(e, vertexWeighted);
    }
    abstract public double getHyperedgeInDegree(int e, boolean vertexWeighted);

    // hyeredge out degrees
    public final List<Double> getHyperedgeOutDegrees() {
        return getHyperedgeOutDegrees(vertexWeighted);
    }
    abstract public List<Double> getHyperedgeOutDegrees(boolean vertexWeighted);
    // hyperedge out degree
    public double getHyperedgeOutDegree(int e) {
        return getHyperedgeOutDegree(e, vertexWeighted);
    }
    abstract public double getHyperedgeOutDegree(int e, boolean vertexWeighted);

    // hyperedge degrees
    public final List<Double> getHyperedgeDegrees() {
        return getHyperedgeDegrees(vertexWeighted);
    }
    abstract public List<Double> getHyperedgeDegrees(boolean vertexWeighted);
    // hyperedge degree
    public double getHyperedgeDegree(int e) {
        return getHyperedgeDegree(e, vertexWeighted);
    }
    abstract public double getHyperedgeDegree(int e, boolean vertexWeighted);



    ///////////////////////////////////////////////////////////////////////////
    //////// methods for calculating vertex degrees
    //////////////////////////////////////////////////////////////////////////
    // tail vertex degrees
    public final List<Double> getTailVertexDegrees() {
        return getTailVertexDegrees(hyperedgeWeighted);
    }
    abstract public List<Double> getTailVertexDegrees(boolean hyperedgeWeighted);
    // tail vertex degree
    public double getTailVertexDegree(int v) {
        return getTailVertexDegree(v, hyperedgeWeighted);
    }
    abstract public double getTailVertexDegree(int v, boolean hyperedgeWeighted);

    // head vertex degrees
    public final List<Double> getHeadVertexDegrees() {
        return getHeadVertexDegrees(hyperedgeWeighted);
    }
    abstract public List<Double> getHeadVertexDegrees(boolean hyperedgeWeighted);
    // head vertex degree
    public double getHeadVertexDegree(int v) {
        return getHeadVertexDegree(v, hyperedgeWeighted);
    }
    abstract public double getHeadVertexDegree(int v, boolean hyperedgeWeighted);

    // vertex degrees
    public List<Double> getVertexDegrees() {
        return getVertexDegrees(hyperedgeWeighted);
    }
    abstract public List<Double> getVertexDegrees(boolean hyperedgeWeighted);
    // vertex degree
    public double getVertexDegree(int v) {
        return getVertexDegree(v, hyperedgeWeighted);
    }
    abstract public double getVertexDegree(int v, boolean hyperedgeWeighted);


    abstract public Matrix getOutAdjacencyMatrix(boolean vertexWeighted);
    abstract public Matrix getInAdjacencyMatrix(boolean vertexWeighted);


    //////////////////////////////////////////////////////////////
    ////////// Methods for adjacency matrix normalization
    //////////////////////////////////////////////////////////////
    public Matrix normalize() {
        return normalize(hyperedgeWeighted, vertexWeighted/*, true*/);
    }
    abstract public Matrix normalize(boolean hyperedgeWeighted, boolean vertexWeighted/*, boolean allowSelfLoop*/);


    public Matrix normalize(boolean H2A) {
        return normalize(H2A, hyperedgeWeighted, vertexWeighted/*, true*/);
    }
    abstract public Matrix normalize(boolean H2A, boolean hyperedgeWeighted, boolean vertexWeighted/*, boolean allowSelfLoop*/);

}
