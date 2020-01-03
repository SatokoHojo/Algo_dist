import io.jbotsim.core.Node;
import io.jbotsim.core.Message;
import io.jbotsim.core.Color;
import java.util.List;
import java.lang.Math;
import java.util.BitSet;

public class MyNode extends Node{
    static int nb_nodes;
    private Node father = null;
    private int color;
    private int father_color;
    private int l;
    private int l_prime;
    private int to_remove = 5;
    private boolean shift = true;

    private void getCorrespColo(int c) {
        switch (c) {
            case 0:
                setColor(Color.BLUE);
                break;
            case 1:
                setColor(Color.RED);
                break;
            case 2:
                setColor(Color.GREEN);
                break;
            case 3:
                setColor(Color.BLACK);
                break;
            case 4:
                setColor(Color.ORANGE);
                break;
            case 5:
                setColor(Color.YELLOW);
                break;
            default:
                throw new IllegalStateException("Unexpected value: " + c);
        }
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

    private int FirstFree(List<Message> l) {
        //here we suppose that there are at most 6 colors
        boolean[] present = {false, false, false, false, false, false};
        for (Message m: l) {
            present[(int)m.getContent()] = true;
        }
        for (int i = 0; i < 6; i++) {
            if (!present[i]) return i;
        }
        throw new IllegalStateException("More than 6 colors among the neighbours !");
    }

    @Override
    public void onStart() {
        /*
         JBotSim executes this method on each node upon initialization
        For now we consider that the father of a node is the neighbour with the smallest id
        */
        color = getID();
        father_color = Math.max(1 - color, 0);
        l = (int) Math.ceil(Math.log(nb_nodes)/Math.log(2));
        l_prime = l+1; //arbitrarily fixed so that l != l_prime
        List<Node> nei = getNeighbors();
        father = this;
        for (Node n: nei) {
            if (father == this || n.getID() > father.getID()) {
                father = n;
            }
        }
    }

    @Override
    public void onSelection() {
        // JBotSim executes this method on a selected node
    }

    @Override
    public void onClock() {
        // JBotSim executes this method on each node in each round
        List<Message> messages = getMailbox();
        if (to_remove >= 3) {
            if (!messages.isEmpty()) {
                //at the first iteration, the mailbox is empty
                for (Message m : messages) {
                    if (m.getSender() == father) {
                        father_color = (int) m.getContent();
                    }
                }
                if (l != l_prime) {
                    color = PosDiff(color, father_color);
                    l_prime = l;
                    l = 1 + (int) Math.ceil(Math.log(l) / Math.log(2));
                } else {
                    getCorrespColo(color);
                    if (shift) {
                        color = father_color;
                        shift = false;
                    } else {
                        if (color == to_remove) {

                        }
                        to_remove = to_remove -1;
                        shift = true;
                    }
                    getCorrespColo(color);
                }
            }
        } else {
            //conflicts
            for (Message m : messages) {
                if ((int) m.getContent() == color && m.getSender().getID() > this.getID()) {
                    color = FirstFree(messages);
                    break;
                }
            }
            getCorrespColo(color);
        }
        sendAll(new Message(color));
    }

    @Override
    public void onMessage(Message message) {
        // JBotSim executes this method on a node every time it receives a message
    }
}