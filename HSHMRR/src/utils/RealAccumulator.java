// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package utils;

public class RealAccumulator implements RealValued {
    protected double value = 0;

    @Override
    public String toString() {
        String str = String.format("%.6f", value);
        return str;
    }

    @Override
    public double getRealValue() {
        // TODO Auto-generated method stub
        return value;
    }

    @Override
    public void setValue(double val) {
        // TODO Auto-generated method stub
        value = val;
    }

    public RealAccumulator (){
        value =0 ;
    }
    public RealAccumulator(double val) {
        value = val;
    }

    public RealAccumulator add(double val){
        value += val;
        return this;
    }

    public void substract(double val) {
        value -= val;
    }

    @Override
    public int compareTo(Object arg0) {
        // TODO Auto-generated method stub
        return (new Double(value)).compareTo(((RealAccumulator)arg0).value);
    }
}
