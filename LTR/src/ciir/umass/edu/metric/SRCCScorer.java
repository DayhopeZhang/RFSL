package ciir.umass.edu.metric;

import ciir.umass.edu.learning.RankList;

public class SRCCScorer extends MetricScorer {

    public SRCCScorer() {this.k = 10;}
    public SRCCScorer(int k) {this.k = k;}
    public SRCCScorer copy() {return new SRCCScorer();}

    public double score(RankList rl) {
        float[] rel = new float[rl.size()];
        for (int i = 0; i < rel.length; ++i)
            rel[i] = rl.get(i).getLabel();

        float[] rank = new float[rl.size()];
        for (int i = 0; i < rel.length; ++i) {
            int count1 = 1;
            int count2 = -1;
            for (int j = 0; j < rel.length; ++j) {
                if (Math.abs(rel[j] - rel[i]) < 1e-10)
                    count2++;
                else if (rel[j] > rel[i])
                    count1++;
            }
            rank[i] = count1 + count2 / 2F;
        }

        double meanX = (1 + rl.size()) / 2.0;
        double meanY = 0;
        for (int i = 0; i < rank.length; ++i)
            meanY += rank[i];
        meanY /= rl.size();

        double numerator = 0, denominator1 = 0, denominator2 = 0;
        for (int i = 0; i < rank.length; ++i) {
            numerator += (i + 1 - meanX) * (rank[i] - meanY);
            denominator1 += (i + 1 - meanX) * (i + 1 - meanX);
            denominator2 += (rank[i] - meanY) * (rank[i] - meanY);
        }

        return numerator / Math.sqrt(denominator1 * denominator2);
    }

    public double[][] swapChange(RankList rl) {
        float[] rel = new float[rl.size()];
        for (int i = 0; i < rel.length; ++i)
            rel[i] = rl.get(i).getLabel();

        float[] rank = new float[rl.size()];
        for (int i = 0; i < rel.length; ++i) {
            int count1 = 1;
            int count2 = -1;
            for (int j = 0; j < rel.length; ++j) {
                if (Math.abs(rel[j] - rel[i]) < 1e-10)
                    count2++;
                else if (rel[j] > rel[i])
                    count1++;
            }
            rank[i] = count1 + count2 / 2F;
        }

        double meanX = (1 + rl.size()) / 2.0;
        double meanY = 0;
        for (int i = 0; i < rank.length; ++i)
            meanY += rank[i];
        meanY /= rl.size();

        double denominator1 = 0, denominator2 = 0;
        for (int i = 0; i < rank.length; ++i) {
            denominator1 += (i + 1 - meanX) * (i + 1 - meanX);
            denominator2 += (rank[i] - meanY) * (rank[i] - meanY);
        }
        double denominator = Math.sqrt(denominator1 * denominator2);

        double[][] change = new double[rl.size()][rl.size()];
        for (int i = 0; i < rl.size(); ++i) {
            for (int j = i + 1; j < rl.size(); ++j) {
                double delta = (i + 1 - meanX) * (rank[j] - meanY) + (j + 1 - meanX) * (rank[i] - meanY)
                        - (i + 1 - meanX) * (rank[i] - meanY) - (j + 1 - meanX) * (rank[j] - meanY);
                delta /= denominator;
                change[i][j] = change[j][i] = delta;
            }
        }

        return change;
    }

    public String name() {return "SRCCScorer";}
}
