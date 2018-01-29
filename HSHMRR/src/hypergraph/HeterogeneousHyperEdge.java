// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class HeterogeneousHyperEdge implements Comparable<HeterogeneousHyperEdge> {
    int heid;
    int nmode;
    // vmode -- vid
    List<List<Integer>> vids;	//

    public HeterogeneousHyperEdge(int heid, int nmode) {
        this(heid, nmode, new ArrayList<List<Integer>>(nmode));
    }

    public HeterogeneousHyperEdge(int heid, int nmode, List<List<Integer>> vids) {
        this.heid = heid;
        this.nmode = nmode;
        this.vids = vids;
        consolidate();
    }

    @Override
    public String toString() {
        String str = "" + heid + "[";
        for (int vmode = 0; vmode < nmode; ++vmode) {
            str += "mode_" + vmode + vids.get(vmode).toString() + ", ";
        }
        str = str.substring(0, str.length()-2) + "]";
        return str;
    }

    public int getNumderOdModes() {
        return nmode;
    }

    public List<List<Integer>> getVertexIDs() {
        return vids;
    }

    public List<Integer> getVertexIDs(int vmode) {
        return vids.get(vmode);
    }

    public int getHyperEdgeID() {
        return heid;
    }

    public void addVertex(int v, int vmode) {
        if (!vids.contains(v))
            vids.get(vmode).add(v);
        consolidate(vmode);
    }

    public void addVertices(List<List<Integer>> vids) {
        for (int vmode = 0; vmode < nmode; ++vmode) {
            addVertices(vids.get(vmode), vmode);
        }
    }

    public void addVertices(List<Integer> vids, int vmode) {
        vids.removeAll( this.vids.get(vmode) );
        this.vids.get(vmode).addAll(vids);
        consolidate(vmode);
    }

    /**
     * Sorting the vids
     */
    public void consolidate() {		// make the sorted
        for (int vmode = 0; vmode < nmode; ++vmode)
            if (!this.vids.isEmpty())	// Dealing with empty hyperedge as just after initialized
                consolidate(vmode);
    }

    public void consolidate(int vmode) {		// make the sorted
        Collections.sort(this.vids.get(vmode));
    }

    public int getHypeEdgeDegree(int vmode) {
        return vids.get(vmode).size();
    }

    @Override
    public boolean equals(Object o) {
        boolean eq = true;
        if (o instanceof HyperEdge) {
//			//Cater to libANN : buildAuthorCollaborationHypergraph; seems to need modification 20150716
//			if (heid != ((HyperEdge) o).heid)
//				eq = false;
//			else {
            if (vids.size() != ((HyperEdge) o).vids.size())
                eq = false;
            else if (nmode !=  ((HeterogeneousHyperEdge) o).nmode ) {
                eq = false;
            } else {
                for (int vmode = 0; vmode < nmode; ++vmode) {
                    if (vids.get(vmode).size() != ((HeterogeneousHyperEdge) o).vids.get(vmode).size() ) {
                        eq = false; break;
                    }
                    if (!this.vids.get(vmode).equals(((HeterogeneousHyperEdge)o).vids.get(vmode))) {
                        eq = false; break;
                    }
                }
            }
//			}
        } else
            eq = false;
        return eq;
    }

    @Override
    public int hashCode() {
        int z = 17;
        int hc = 1;
        for (int vmode = 0; vmode < nmode; ++vmode) {
            for (int i = 0; i < vids.get(vmode).size(); ++i) {
                hc = hc * z + vids.get(vmode).get(i);
            }
        }
        return hc;
    }

    @Override
    public int compareTo(HeterogeneousHyperEdge arg0) {
        // TODO Auto-generated method stub
        return ((Integer) heid).compareTo(arg0.heid);
    }
}
