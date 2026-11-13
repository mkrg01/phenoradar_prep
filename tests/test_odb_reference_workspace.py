"""Reference lists are writable only in each mapping job's private workspace."""
import pytest
from run_odb_chunk import stage_reference_data


def test_todo_sort_and_new_lists_are_private_and_resumable(tmp_path):
    source = tmp_path / 'source'
    todo = source / 'eukaryota/todo/g3193_n8ssc.todo'
    todo.parent.mkdir(parents=True)
    todo.write_text('beta\nalpha\n')
    database = source / 'eukaryota/sequences.db'
    database.write_bytes(b'large immutable database')
    for name in ('first', 'second'):
        dest = stage_reference_data(source, tmp_path / name)
        local = dest / todo.relative_to(source)
        assert not local.samefile(todo)
        # Reproduce Orthologer's in-place sort, including repeated setup on retry.
        local.write_text(''.join(sorted(local.read_text().splitlines(keepends=True))))
        (local.parent / 'generated.todo').write_text('new\n')
        stage_reference_data(source, dest)
        assert local.read_text() == 'alpha\nbeta\n'
        assert (dest / database.relative_to(source)).is_symlink()
    assert todo.read_text() == 'beta\nalpha\n'
    assert not (todo.parent / 'generated.todo').exists()


@pytest.mark.parametrize('link', ['symlink', 'hardlink'])
def test_existing_shared_todo_is_rejected(tmp_path, link):
    source = tmp_path / 'source'; source.mkdir()
    todo = source / 'node.todo'; todo.write_text('beta\nalpha\n')
    dest = tmp_path / 'dest'; dest.mkdir()
    target = dest / todo.name
    if link == 'symlink': target.symlink_to(todo)
    else: target.hardlink_to(todo)
    with pytest.raises(ValueError, match='shares the frozen reference'):
        stage_reference_data(source, dest)
    assert todo.read_text() == 'beta\nalpha\n'
