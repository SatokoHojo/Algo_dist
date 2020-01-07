import io.jbotsim.core.Color;
import io.jbotsim.core.Node;
import io.jbotsim.core.Message;
import java.util.ArrayList;
import java.lang.Math;
import java.util.BitSet;
import java.util.List;

public class GeneralNode extends MyNode {
    private ArrayList<Node> fathers;
    private int[] colors;
    private int[] father_colors;
    private boolean shift = true;
    private int final_color = -1;
    private boolean reduction_started = false;
    private boolean printed_res = false;

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
        return (max + 1);
    }

    void sendColors(){
        int[] colors_to_send = new int[colors.length];
        for (int i = 0 ; i < colors.length; i++){
            colors_to_send[i] = colors[i];
        }
        sendAll(new Message(colors_to_send));
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


    private void color6Dloop(){
        for (int i = 0; i < delta; i++) {
            //System.out.print("dans color father_colors");
            //System.out.println("id : " + getID() + ", color : " + colors[i] + ", father's : " + father_colors[i] + ", father : " + fathers.get(i).getID());
            colors[i] = PosDiff(colors[i], father_colors[i]);
            l_prime = l;
            l = 1 + (int) Math.ceil(log2(l));
        }
        sendColors();
    }

    private void reduce3to6D(List<Message> messages){
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


    private void reduceColors(List<Message> messages){
        if (to_remove >= delta+1) {
            if (final_color == to_remove) {
                final_color = FirstFree(messages);
            }
            to_remove = to_remove -1;
        }
        sendAll(new Message(final_color) );
    }

    @Override
    public void onClock() {
        if (! finished){

            //partie 2 : reduction de 3D à D+1 couleurs
            if (reduction_started) {
                List<Message> messages = getMailbox();
                reduceColors(messages);
                if (to_remove == delta){
                    getCorrespColor(final_color);
                    finished = true;
                    System.out.print("final : id : " + getID() + ", final color : " + final_color + "\n");
                }
            }

            //Partie 1 : 3D coloration
            else {
                List<Message> messages = getMailbox();
                for (Message m : messages) {
                    if (m.getSender().getID() < getID()) {
                        int index = fathers.indexOf(m.getSender());
                        father_colors[index] = ((int[]) m.getContent())[index];
                    }
                }
                if (l != l_prime) {
                    color6Dloop();
                } else if (to_remove >= 3) { //passer de 6 à 3 couleur par arbres
                    reduce3to6D(messages);
                } else { // on a fini notre 3D coloration, on passe à la 2eme partie
                    int tmp = 0;
                    for (int i = 0; i < delta; i++) tmp += (int) Math.pow(3, i) * colors[i];
                    final_color = tmp;
                    to_remove = (int) Math.pow(3, delta) - 1;
                    reduction_started = true;
                    System.out.println("id : " + getID() + ", color : " + final_color + " nb colors : " + to_remove);
                    sendAll(new Message(final_color) );
                }
            }
        }

        /*


        // JBotSim executes this method on each node in each round
        List<Message> messages = getMailbox();
        if (to_remove >= 3 && !reduction_started) { //1ere partie
            System.out.print("id : " + getID() + ", current colors : " + colors[0] + ", " + colors[1] + "\n");
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
                        System.out.println("id : " + getID() + ", color : " + colors[i] + ", father's : " + father_colors[i] + ", father : " + fathers.get(i).getID());
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
            System.out.print("id : " + getID() + ", new colors : " + colors[0] + ", " + colors[1] + "\n");
            sendAll(new Message(colors));



        } else { // reduce palettte
            // a priori no conflicts
            if (!reduction_started) {
                int tmp = 0;
                for (int i = 0; i < delta; i++) tmp += (int) Math.pow(3, i)*colors[i];
                final_color = tmp;
                System.out.println("id : " + getID() + ", color : " + final_color);
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
                System.out.print("id : " + getID() + ", final color : " + final_color + "\n");
                printed_res = true;
            }
            sendAll(new Message(final_color));
        }*/
    }
}
