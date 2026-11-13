from docutils import nodes
from sphinx.util.docutils import SphinxDirective


class Collapsible(nodes.General, nodes.Element):
    pass


class CollapsedSource(SphinxDirective):
    has_content = True
    required_arguments = 1
    final_argument_whitespace = True
    option_spec = {}

    def run(self):
        wrapper = Collapsible()
        wrapper["title"] = self.arguments[0]

        self.state.nested_parse(
            self.content,
            self.content_offset,
            wrapper,
        )
        return [wrapper]

def visit_collapsible_html(self, node):
    self.body.append(
        '<details class="collapsed-source">'
        '<summary style="'
        'cursor:pointer;'
        'font-size:1.15em;'
        'font-weight:700;'
        'color:#174a6e;'
        'background:#eaf3f8;'
        'border:1px solid #9bbccc;'
        'border-radius:4px;'
        'padding:0.6em 0.8em;'
        '">'
        f'{self.encode(node["title"])}'
        '</summary>\n'
    )


def visit_collapsible_other(self, node):
    # For non-HTML builders, render the contents without the wrapper.
    pass

def depart_collapsible_html(self, node):
    self.body.append("</details>\n")

def depart_collapsible_other(self, node):
    pass


def setup(app):
    app.add_node(
        Collapsible,
        html=(visit_collapsible_html, depart_collapsible_html),
        latex=(visit_collapsible_other, depart_collapsible_other),
        text=(visit_collapsible_other, depart_collapsible_other),
        man=(visit_collapsible_other, depart_collapsible_other),
        texinfo=(visit_collapsible_other, depart_collapsible_other),
    )
    app.add_directive("collapsed-source", CollapsedSource)

    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
