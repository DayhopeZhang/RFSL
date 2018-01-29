// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class HyperEdge implements Comparable<HyperEdge> {

    int heid;
    List<Integer> vids;	//

    public HyperEdge(int heid) {
        this(heid, new ArrayList<Integer>());
    }

    public HyperEdge(int heid, List<Integer> vids) {
        this.heid = heid;
        this.vids = vids;
        consolidate();
    }

    public List<Integer> getVertexIDs() {
        return vids;
    }

    public int getHyperEdgeID() {
        return heid;
    }

    public void addVertex(int v) {
        if (!vids.contains(v))
            vids.add(v);
        consolidate();
    }

    @Override
    public String toString() {
        String str = "";
        str += heid + vids.toString();
        return str;
    }

//	public void addVertexAndConsolidate(int v) {
//		addVertex(v);
//		consolidate();
//	}

    public void addVertices(List<Integer> vids) {
        for (int v : vids) {
            if (this.vids.contains(v))
                this.vids.add(v);
        }
        consolidate();
    }

//	public void addVerticesAndConsolidate(List<Integer> vids) {
//		addVertices(vids);
//		consolidate();
//	}

    /**
     * Sorting the vids
     */
    public void consolidate() {		// make the sorted
        Collections.sort(vids);
    }

    public int getHypeEdgeDegree() {
        return vids.size();
    }

    public boolean equals(Object o) {
        boolean eq = true;
        if (o instanceof HyperEdge) {
//			//Cater to libANN : buildAuthorCollaborationHypergraph; seems to need modification 20150716
//			if (heid != ((HyperEdge) o).heid)
//				eq = false;
//			else {
            if (vids.size() != ((HyperEdge) o).vids.size())
                eq = false;
            else {
                if (vids.equals(((HyperEdge) o).vids))
                    eq = false;
//					for (int i = 0; i < vids.size(); ++i) {
//						if (vids.get(i) != ((HyperEdge) o).vids.get(i)) {
//							eq = false; break;
//						}
//					}
            }
//			}
        } else
            eq = false;
        return eq;
    }

    public int hashCode() {
        int z = 17;
        int hc = 1;
        for (int i = 0; i < vids.size(); ++i) {
            hc = hc * z + vids.get(i);
        }
        return hc;
    }

    @Override
    public int compareTo(HyperEdge arg0) {
        // TODO Auto-generated method stub
        return ((Integer) heid).compareTo(arg0.heid);
    }
}
