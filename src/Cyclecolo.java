import io.jbotsim.core.Topology;
import io.jbotsim.ui.JViewer;
import io.jbotsim.core.Node;

public class Cyclecolo {
    public static void main(String[] args){
        MyTopology tp = new MyTopology();
        tp.setDefaultNodeModel(GeneralNode.class);
        new JViewer(tp);
        //tp.start();
        //System.out.println(tp.getNodes().size());
    }
}