from sqlalchemy import create_engine, MetaData

echo = False

dbfile = '/home/nhoffman/src/Taxonomy/test_output/taxonomy_test.db'
engine = create_engine('sqlite:///%s' % dbfile, echo=echo)

meta = MetaData()
meta.bind = engine
meta.reflect()

for t in meta.sorted_tables:
    print t.name

nodes = meta.tables['nodes']
names = meta.tables['names']



    # s = nodes.select()
    # print s
    # result = s.execute()
    # for row in result:
    #     print row
    #     break
