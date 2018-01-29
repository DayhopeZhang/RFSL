// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package utils;

public interface RealValued extends Comparable {
    public double getRealValue();

    public void setValue(double val);
}
