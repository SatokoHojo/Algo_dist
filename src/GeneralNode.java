import io.jbotsim.core.Node;
import io.jbotsim.core.Message;
import java.util.ArrayList;
import java.lang.Math;
import java.util.List;

public class GeneralNode extends MyNode {
    private ArrayList<Node> fathers;
    private int[] colors;
    private int[] father_colors;
    private boolean shift;
    private boolean reduction_started;

    //this function is called only by the topology before the beginning of the clock
    //on all nodes in a sequatial way.
    public void init_delta() {
        int degree = getNeighbors().size();
        if (degree > delta){
            delta = degree;
        }
    }


    @Override
    public void onStart() {

        super.onStart();
        shift = true;
        reduction_started = false;

        colors = new int[delta];
        for (int i = 0; i < delta; i++) colors[i] = getID();

        // initializing the father
        //For now we consider that any neighbour with a smaller id is a father
        father_colors = new int[delta];
        fathers = new ArrayList<>(delta);

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
        sendColors();
    }



    /* A OPTIMISER*/
    //here l is the list of message containing the color (table) of each neighbor
    //index is the number of the 1-orientation we ar focusing on
    private int FirstFree(List<Message> l, int index) {
        // careful with the number of neighbours, or this function will explode
        int max = 0;
        for (Message m : l) {
            max = Math.max(max, ((int[]) m.getContent())[index]);
        }
        boolean[] present = new boolean[max+1];
        for (Message m: l) {
            present[((int[])m.getContent())[index]] = true;
        }
        for (int i = 0; i < max+1; i++) {
            if (!present[i]) return i;
        }
        return (max + 1);
    }


    //here l is the list of message containing the color of each neighbor
    protected int FirstFree(List<Message> l) {
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
        return (max + 1);
    }
    /**/

    void sendColors(){
        int[] colors_to_send = new int[colors.length];
        System.arraycopy(colors, 0, colors_to_send, 0, colors.length);
        sendAll(new Message(colors_to_send));
    }

    //main part of the algorithm that color with 6 colors each 1-orientation
    private void color6Dloop(){
        for (int i = 0; i < delta; i++) {
            colors[i] = PosDiff(colors[i], father_colors[i]);
        }
        l_prime = l;
        l = 1 + (int) Math.ceil(log2(l));
        sendColors();
    }

    //main part of the algorithm reduce from 6 to 3 the number of colors for each 1-orientation
    private void reduce3to6D(List<Message> messages){ //messages is the list of the colors (table) of the neighbors
        if (shift) {
            for (int i = 0; i < delta; i++) {
                colors[i] = father_colors[i];
                if (fathers.get(i).getID() == getID()){
                    father_colors[i] = Math.max(1-colors[i], 0);
                }
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
            shift = true;
        }
        sendColors();
    }


    @Override
    public void onClock() {
        if (! finished){
            List<Message> messages = getMailbox(); //messages is the list of the colors of the neighbors

            //PART 2 : reducePalette from 3^delta colors to delta+1
            if (reduction_started) {
                reducePalette(messages);

                //End
                if (to_remove == delta){
                    getCorrespColor(myColor);
                    System.out.print("Result : node " + getID() + ", colored " + myColor + "\n");
                    finished = true;
                }
            }

            //PART 1 : 3^delta coloration (3-coloration for each 1-coloration)
            else {
                //getting the colors (tables) of the neighbors
                for (Message m : messages) {
                    if (m.getSender().getID() < getID()) {
                        int index = fathers.indexOf(m.getSender());
                        father_colors[index] = ((int[]) m.getContent())[index];
                    }
                }
                //PART 1.A 6-coloration for each 1-coloration
                if (l != l_prime) {
                    color6Dloop();

                //PART 1.B Reduce from 6 to 3 colors for each 1-coloration
                } else if (to_remove >= 3) {
                    reduce3to6D(messages);

                } else {
                    // Here the graphes is 3^delta-colored :
                    // preparation for reducepalette from 3^delta to delta +1
                    for (int i = 0; i < delta; i++) myColor += (int) Math.pow(3, i) * colors[i];
                    to_remove = (int) Math.pow(3, delta) - 1;
                    //System.out.println("id : " + getID() + ", color : " + myColor + " nb colors : " + to_remove);
                    sendAll(new Message(myColor) );
                    reduction_started = true;
                }
            }
        }
    }
}
