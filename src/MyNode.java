import io.jbotsim.core.Color;
import io.jbotsim.core.Node;

import java.util.BitSet;

public abstract class MyNode extends Node {

    static int nb_nodes = -1;
    static int delta = -1;
    protected int l = -1;
    protected int l_prime;
    protected int to_remove = 5;
    protected  boolean finished = false;

    public void init() {
        nb_nodes++;
        int degree = getNeighbors().size();
        if (degree > delta){
            delta = degree;
        }
    }

    protected double log2(int x) {
        return Math.log(x)/Math.log(2);
    }

    protected BitSet binary(int n, int size) {
        BitSet ret = new BitSet(size);
        int index = 0;
        int tmp = n;
        while (tmp != 0) {
            if (tmp%2 == 1) ret.set(index);
            index = index+1;
            tmp = tmp/2;
        }
        return ret;
    }

    protected int PosDiff(int c, int cf) {
        int m = Math.max(c,cf);
        int max_pow = 0;
        int tmp = 2;
        while (m > tmp-1) {
            max_pow = max_pow+1;
            tmp = tmp*2;
        }
        int size = max_pow +1;
        BitSet bin_c = binary(c, size);
        BitSet bin_cf = binary(cf, size);
        BitSet diff = (BitSet) bin_c.clone();
        diff.xor(bin_cf);
        if (diff.isEmpty()) throw new IllegalStateException("Same colors !");
        int p = 0;
        while (!diff.get(p)) { //while the bits are identical
            p++;
        }
        int bin_p = bin_c.get(p) ? 1 : 0;
        return 2*p + bin_p;
    }

    protected void getCorrespColor(int c) {
        int max_color = delta;
        if (c > max_color) {
            throw new IllegalStateException("Unexpected value: " + c);
        } else {
            if(max_color == 0) setColor(new Color(0));
            else {
                int col = (int) Math.floor(c * (255. / max_color));
                //System.out.print(col + " ");
                setColor(new Color(col, 255 - col, 0));
            }
        }
    }
}
