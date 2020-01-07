import io.jbotsim.ui.JViewer;

public class CycleColo {
    public static void main(String[] args){
        CycleTopology tp = new CycleTopology();
        tp.setDefaultNodeModel(CycleNode.class);
        new JViewer(tp);
    }
}