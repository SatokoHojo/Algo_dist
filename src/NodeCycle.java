import io.jbotsim.core.Node;
import io.jbotsim.core.Message;
import io.jbotsim.core.Color;
import java.util.List;
import java.lang.Math;
import java.util.BitSet;

public class NodeCycle extends MyNode {
    private Node father = null;
    private int color;
    private int father_color;
    private int succ_color;


/*
    protected void getCorrespColor(int c) {
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
*/

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
        nb_nodes = Math.max(nb_nodes, color);
        father = this;
    }

    @Override
    public void onSelection() {
        // JBotSim executes this method on a selected node
    }

    @Override
    public void onClock() {
        // JBotSim executes this method on each node in each round
        if (l < 0) {
            l = (int) Math.ceil(log2(nb_nodes));
            l_prime = l+1; //arbitrarily fixed so that l != l_prime
            // initializing the father
            List<Node> nei = getNeighbors();
            for (Node n: nei) {
                if (father == this || n.getID() > father.getID()) {
                    father = n;
                }
            }
        }
        List<Message> messages = getMailbox();
        if (to_remove >= 3) {
            if (!messages.isEmpty()) {
                //at the first iteration, the mailbox is empty
                for (Message m : messages) {
                    if (m.getSender() == father) father_color = (int) m.getContent();
                    else succ_color = (int) m.getContent();
                }
                if (l != l_prime) {
                    color = PosDiff(PosDiff(father_color, color),PosDiff(color, succ_color));
                    l_prime = l;
                    l = 1 + (int) Math.ceil(log2(1 + (int) Math.ceil(log2(l))));
                } else {
                    getCorrespColor(color);
                    if (color == to_remove) {
                        //ReducePalette
                        color = FirstFree(messages);
                    }
                    to_remove = to_remove -1;
                    getCorrespColor(color);
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
            getCorrespColor(color);
        }
        sendAll(new Message(color));
    }

    @Override
    public void onMessage(Message message) {
        // JBotSim executes this method on a node every time it receives a message
    }
}
