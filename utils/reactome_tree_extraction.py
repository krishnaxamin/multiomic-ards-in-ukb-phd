""" Get all ancestors for all Reactome annotations. """
from pandas import read_csv, DataFrame
from collections import defaultdict

hierarchy = read_csv('~/data/external/reactome/ReactomePathwaysRelation.txt', sep='\t', header=None).set_axis(['parent', 'child'], axis=1)
hierarchy_human = hierarchy[hierarchy['child'].str.contains('HSA')].copy()
child_to_parent_map = dict(zip(hierarchy_human['child'], hierarchy_human['parent']))


def get_all_ancestors(child, ancestors=None):

    if ancestors is None:
        ancestors = {child}

    if child not in child_to_parent_map.keys():
        return set()

    parent = child_to_parent_map[child]
    ancestors.add(parent)

    return ancestors | get_all_ancestors(parent, ancestors=ancestors)


id_ancestors = defaultdict(list)
for child in hierarchy_human['child'].unique():
    id_ancestors[child] = list(get_all_ancestors(child))

id_ancestors_df = DataFrame({'id': list(id_ancestors.keys()),
                             'ancestors': [x for x in id_ancestors.values()]})
id_ancestors_df = id_ancestors_df.explode('ancestors', ignore_index=True)
id_ancestors_df.to_csv('~/data/external/reactome/reactome/reactome_id_ancestor_relations.csv', index=False)
