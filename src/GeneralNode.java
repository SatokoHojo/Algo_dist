import io.jbotsim.core.Node;
import io.jbotsim.core.Message;
import java.util.ArrayList;
import java.lang.Math;
import java.util.BitSet;
import java.util.List;

public class GeneralNode extends Node {
    static int nb_nodes = 10; // assuming for now we calculate it somehow before starting
    static int delta = 2; // same
    private ArrayList<Node> fathers;
    private int[] colors;
    private int[] father_colors;
    private int l = -1;
    private int l_prime;
    private boolean init = true;
    private int to_remove = 5;
    private boolean shift = true;
    private int final_color = -1;
    private boolean reduction_started = false;
    private boolean printed_res = false;

    private double log2(int x) {
        return Math.log(x)/Math.log(2);
    }

    private BitSet binary(int n, int size) {
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

    private int PosDiff(int c, int cf) {
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

    private int FirstFree(List<Message> l, int index) {
        // careful with the number of neighbours, or this function will explode
        int max = 0;
        for (Message m : l) {
            max = Math.max(max, ((int[]) m.getContent())[index]);
        }
        boolean[] present = new boolean[max];
        for (Message m: l) {
            present[((int[])m.getContent())[index]] = true;
        }
        for (int i = 0; i < max+1; i++) {
            if (!present[i]) return i;
        }
        throw new IllegalStateException("More than " + max + " colors among the neighbours !");
    }

    private int FirstFree(List<Message> l) {
        // careful with the number of neighbours, or this function will explode
        int max = 0;
        for (Message m : l) {
            max = Math.max(max, (int) m.getContent());
        }
        boolean[] present = new boolean[max+1];
        for (Message m: l) {
            present[(int)m.getContent()] = true;
        }
        for (int i = 0; i < max+1; i++) {
            if (!present[i]) return i;
        }
        throw new IllegalStateException("More than " + max + " colors among the neighbours !");
    }

    @Override
    public void onStart() {
        /*
         JBotSim executes this method on each node upon initialization
        For now we consider that any neighbour with a smaller id is a father
        */
        colors = new int[delta];
        for (int i = 0; i < delta; i++) colors[i] = getID();

        father_colors = new int[delta];

        fathers = new ArrayList<Node>(delta);

        l = (int) Math.ceil(log2(nb_nodes));
        l_prime = l+1; //arbitrarily fixed so that l != l_prime
    }

    @Override
    public void onClock() {
        // JBotSim executes this method on each node in each round
        if (init) { // because the list of neighbours isn't fixed in onStart
            List<Node> nei = getNeighbors();
            int current_nb = 0;
            for (Node n : nei) {
                if (n.getID() < getID()) {
                    fathers.add(n);
                    current_nb++;
                }
            }
            for (int i = current_nb; i < delta; i++) {
                fathers.add(this);
                father_colors[i] = Math.max(1-getID(), 0);
            }
            init = false;
        }
        List<Message> messages = getMailbox();
        if (to_remove >= 3 && !reduction_started) {
            if (!messages.isEmpty()) {
                //at the first iteration, the mailbox is empty
                for (Message m : messages) {
                    if (m.getSender().getID() < getID()) {
                        int index = fathers.indexOf(m.getSender());
                        father_colors[index] = ((int[]) m.getContent())[index];
                    }
                }
                if (l != l_prime) {
                    for (int i = 0; i < delta; i++) {
                        colors[i] = PosDiff(colors[i], father_colors[i]);
                        l_prime = l;
                        l = 1 + (int) Math.ceil(log2(l));
                    }
                } else {
                    if (shift) {
                        for (int i = 0; i < delta; i++) {
                            colors[i] = father_colors[i]; //should we shift only if the father is not the node itself ?
                        }
                        shift = false;
                    } else {
                        for (int i = 0; i < delta; i++) {
                            if (colors[i] == to_remove) {
                                //ReducePalette
                                colors[i] = FirstFree(messages, i);
                            }
                        }
                        to_remove = to_remove - 1;
                    }
                }
            }
            sendAll(new Message(colors));
        } else {
            // a priori no conflicts
            if (!reduction_started) {
                int tmp = 0;
                for (int i = 0; i < delta; i++) tmp += (int) Math.pow(3, i)*colors[i];
                final_color = tmp;
                to_remove = (int) Math.pow(3, delta) -1;
                reduction_started = true;
            }
            if (to_remove >= delta) {
                if (final_color == to_remove) {
                    final_color = FirstFree(messages);
                }
                to_remove = to_remove -1;
            } else if (!printed_res) {
                //prints once the final color
                System.out.print("id : " + getID() + ", color : " + final_color + "\n");
                printed_res = true;
            }
            sendAll(new Message(final_color));
        }
    }
}
