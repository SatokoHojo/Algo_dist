import io.jbotsim.core.Color;
import io.jbotsim.core.Message;
import io.jbotsim.core.Node;
import java.util.BitSet;
import java.util.List;

public abstract class MyNode extends Node {

    public static int nb_nodes = 0;
    public static int delta = 0;
    protected int l;
    protected int l_prime;
    protected int to_remove;
    protected  boolean finished;
    protected int myColor;


    @Override
    public void onStart() {
        myColor = 0;
        to_remove = 5;
        finished = false;
        l = (int) Math.ceil(log2(nb_nodes));
        l_prime = l+1;//arbitrarily fixed so that l != l_prime
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


    protected abstract int FirstFree(List<Message> l) ;

    //remove all the colors bigger than delta+1 (leaving delta+1 colors)
    //for the cycle, that corresponds to the color 5,4,3 (leaving the 3 colors 0,1,2)
    protected void reducePalette(List<Message> messages){ //ici les messsages sont des couleurs pour chaque voisins
        if (to_remove > delta) {
            if (myColor == to_remove) {
                myColor = FirstFree(messages);
            }
            to_remove = to_remove -1;
        }
        sendAll(new Message(myColor) );
    }

    protected void getCorrespColor(int c) {
        int max_color = delta;
        if (c > max_color) {
            throw new IllegalStateException("Unexpected value: " + c);
        } else {
            if(max_color == 0) setColor(new Color(0));
            else {
                int col = (int) Math.floor(c * (255. / max_color));
                setColor(new Color(col, 255 - col, 0));
            }
        }
    }
}
