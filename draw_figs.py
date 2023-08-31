import plotly.express as px
from pandas import DataFrame


def draw_fig(intervals: list, height=200, width=1000, to_html=True):
    # intervals: [{"contig": contig, "start": start, "end": end, "type": "Phage"}]
    if not intervals:  # empty list
        return ''
    data = DataFrame(intervals)
    data['length'] = data['end'] - data['start'] + 1
    data['color'] = [1 if x == 'Phage' else 0 for x in data['type']]
    #   contig  start  end       type  length
    # 0      1      1   20  Bacterial      20
    # 1      1     21   30      Phage      10
    # 2      1     31   40  Bacterial      10
    fig = px.bar(data, x='length', y='contig', color='color',
                 hover_data=['start', 'end', 'type'], orientation='h', height=height, width=width)
    fig.update_layout(barmode='stack',
                      coloraxis={'colorbar': {'title': {'text': 'type'},
                                              'tickmode': 'array',
                                              'tickvals': [0, 1],
                                              'ticktext': ['Bacterial', 'Phage']},
                                 'colorscale': [[0, '#0086F0'], [1, '#00F010']]},
                      yaxis={'categoryorder': 'category ascending'})

    # return a div containing the fig
    if to_html:
        return fig.to_html(include_plotlyjs=False, full_html=False)
    else:
        return fig


if __name__ == '__main__':
    d = [{"contig": "1", "start": 1, "end": 20, "type": "Bacterial"},
         {"contig": "1", "start": 21, "end": 30, "type": "Phage"},
         {"contig": "1", "start": 31, "end": 40, "type": "Bacterial"}]
    draw_fig(d, to_html=False).show()
