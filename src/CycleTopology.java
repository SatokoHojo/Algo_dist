import io.jbotsim.core.Node;
import io.jbotsim.core.Topology;

import java.util.List;

public class CycleTopology extends Topology {

    @Override
    public void start() {
        List<Node> nodeList = getNodes();
        MyNode.nb_nodes = nodeList.size();
        MyNode.delta = 2;
        super.start();
    }
}
