import io.jbotsim.core.Node;
import io.jbotsim.core.Message;
import java.util.List;
import java.lang.Math;

public class CycleNode extends MyNode {
    private Node father = null;
    private int father_color;
    private int succ_color;



    protected int FirstFree(List<Message> l) {
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
        if(nb_nodes < 3){
            myColor = getID();
            getCorrespColor(myColor);
            System.out.print("Result : node " + getID() + ", colored " + myColor + "\n");
            finished = true;
        } else {
            super.onStart();
            // initializing the father
            //For now we consider that the father of a node is the neighbour with the smallest id
            List<Node> nei = getNeighbors();
            for (Node n : nei) {
                if (father == null || n.getID() > father.getID()) {
                    father = n;
                }
            }
            myColor = getID();
            sendAll(new Message(myColor));
        }
    }


    @Override
    public void onClock() {

        if (! finished) {
            List<Message> messages = getMailbox();

            //PART 1 : 6-coloration of a cycle
            for (Message m : messages) {
                if (m.getSender() == father) father_color = (int) m.getContent();
                else succ_color = (int) m.getContent();
            }
            if (l != l_prime) {
                myColor = PosDiff(PosDiff(father_color, myColor),PosDiff(myColor, succ_color));
                l_prime = l;
                l = 1 + (int) Math.ceil(log2(1 + (int) Math.ceil(log2(l))));
            }

            //PART 2 : reducePalette from 6 colors to 3
            else{
                reducePalette(messages);

                //End
                if (to_remove == delta) {
                    getCorrespColor(myColor);
                    System.out.print("Result : node " + getID() + ", colored " + myColor + "\n");
                    finished = true;
                }
            }
            sendAll(new Message(myColor));
        }
    }
}
