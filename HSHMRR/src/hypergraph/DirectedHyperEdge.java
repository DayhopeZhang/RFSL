// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package hypergraph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DirectedHyperEdge implements Comparable<DirectedHyperEdge> {
    int heid;
    List<Integer> tail_vids;	//
    List<Integer> head_vids;	//

    public DirectedHyperEdge(int heid) {
        this(heid, new ArrayList<Integer>(), new ArrayList<Integer>());
    }

    public DirectedHyperEdge(int heid, List<Integer> tail_vids, List<Integer> head_vids) {
        this.heid = heid;
        this.tail_vids = tail_vids;
        this.head_vids = head_vids;
        consolidate();
    }

    @Override
    public String toString() {
        String str = "" + heid + "[tail" + tail_vids.toString() + ", head" + head_vids.toString() + "]";
        return str;
    }

    public final List<List<Integer>> getVertexIDs() {
        List<List<Integer>> vids = new ArrayList<List<Integer>>(2);
        vids.add(tail_vids);
        vids.add(head_vids);
        return vids;
    }

    public final List<Integer> getTailVertexIDs() {
        return tail_vids;
    }

    public final List<Integer> getHeadVertexIDs() {
        return head_vids;
    }

    public int getHyperEdgeID() {
        return heid;
    }

    public void addTailVertex(int v) {
        if (!tail_vids.contains(v))
            tail_vids.add(v);
        consolidate();
    }

    public void addHeadVertex(int v) {
        if (!head_vids.contains(v))
            head_vids.add(v);
        consolidate();
    }

//	public void addVertexAndConsolidate(int v) {
//		addVertex(v);
//		consolidate();
//	}

    public void addVertices(List<Integer> tail_vids, List<Integer> head_vids) {
        tail_vids.removeAll(this.tail_vids);
        this.tail_vids.addAll(tail_vids);
        head_vids.removeAll(this.head_vids);
        this.head_vids.addAll(head_vids);
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
        Collections.sort(tail_vids);
        Collections.sort(head_vids);
    }

    public int getHypeEdgeOutDegree() {
        return tail_vids.size();
    }

    public int getHyperEdgeInDegree() {
        return head_vids.size();
    }

    public boolean equals(Object o) {
        boolean eq = true;
        if (o instanceof DirectedHyperEdge) {
            if (tail_vids.size() != ((DirectedHyperEdge) o).tail_vids.size())
                eq = false;
            else {
                if (tail_vids.equals(((DirectedHyperEdge) o).tail_vids))
                    eq = false;
            }
            if (head_vids.size() != ((DirectedHyperEdge) o).head_vids.size())
                eq = false;
            else {
                if (head_vids.equals(((DirectedHyperEdge) o).head_vids))
                    eq = false;
            }
        } else
            eq = false;
        return eq;
    }

    public int hashCode() {
        int z = 17;
        int hc = 1;
        for (int i = 0; i < tail_vids.size(); ++i) {
            hc = hc * z + tail_vids.get(i);
        }
        for (int i = 0; i < head_vids.size(); ++i) {
            hc = hc * z + head_vids.get(i);
        }
        return hc;
    }

    @Override
    public int compareTo(DirectedHyperEdge arg0) {
        // TODO Auto-generated method stub
        return ((Integer) heid).compareTo(arg0.heid);
    }
}
