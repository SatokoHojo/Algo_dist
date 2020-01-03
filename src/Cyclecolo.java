import io.jbotsim.core.Topology;
import io.jbotsim.ui.JViewer;
import io.jbotsim.core.Node;

public class Cyclecolo {
    public static void main(String[] args){
        Topology tp = new Topology();
        tp.setDefaultNodeModel(NodeCycle.class);
        new JViewer(tp);
        tp.start();
        //System.out.println(tp.getNodes().size());
    }
}