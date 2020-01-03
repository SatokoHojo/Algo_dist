import io.jbotsim.core.Topology;
import io.jbotsim.ui.JViewer;
import io.jbotsim.core.Node;

public class Cyclecolo {
    public static void main(String[] args){
        Topology tp = new Topology();
        tp.setDefaultNodeModel(MyNode.class);
        new JViewer(tp);
        tp.start();
    }
}