import io.jbotsim.core.Topology;
import io.jbotsim.core.Node;
import io.jbotsim.ui.JViewer;

import java.util.List;

public class MyTopology extends Topology {

    public MyTopology () {
        setDefaultNodeModel(MyNode.class);
    }

    @Override
    public void start() {
        List<Node> nodeList = getNodes();
        for (Node node: nodeList ) {
            MyNode myNode = (MyNode) node;
            myNode.init();
        }
        super.start();
    }
}
