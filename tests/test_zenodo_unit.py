from __future__ import annotations

from pathlib import Path
from typing import Any, Iterator
from unittest.mock import MagicMock, patch

import pytest

from vcfcache.integrations import zenodo


class _FakeStreamResponse:
    def __init__(self, chunks: list[bytes], status_code: int = 200):
        self._chunks = chunks
        self.status_code = status_code

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size: int = 1 << 20) -> Iterator[bytes]:
        yield from self._chunks

    def __enter__(self) -> "_FakeStreamResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None


def test_download_doi_downloads_first_file_and_honors_sandbox(tmp_path: Path):
    record = {
        "files": [
            {
                "links": {
                    "self": "https://sandbox.zenodo.org/api/files/bucket/test.tar.gz"
                }
            }
        ]
    }

    record_resp = MagicMock()
    record_resp.json.return_value = record
    record_resp.raise_for_status.return_value = None

    stream_resp = _FakeStreamResponse([b"abc", b"def"])

    with patch("vcfcache.integrations.zenodo.requests.get") as get_mock:
        get_mock.side_effect = [record_resp, stream_resp]
        dest = tmp_path / "cache.tar.gz"
        zenodo.download_doi("10.5072/zenodo.12345", dest, sandbox=True)

        assert dest.read_bytes() == b"abcdef"
        first_call_url = get_mock.call_args_list[0].args[0]
        assert first_call_url.startswith(zenodo.ZENODO_SANDBOX_API)


def test_download_doi_raises_if_no_files(tmp_path: Path):
    record_resp = MagicMock()
    record_resp.json.return_value = {"files": []}
    record_resp.raise_for_status.return_value = None

    with patch("vcfcache.integrations.zenodo.requests.get", return_value=record_resp):
        with pytest.raises(zenodo.ZenodoError):
            zenodo.download_doi("10.5281/zenodo.12345", tmp_path / "x.tar.gz")


def test_create_deposit_posts_to_correct_endpoint():
    dep_resp = MagicMock()
    dep_resp.json.return_value = {"id": 1}
    dep_resp.raise_for_status.return_value = None

    with patch("vcfcache.integrations.zenodo.requests.post", return_value=dep_resp) as post_mock:
        dep = zenodo.create_deposit("tok", sandbox=True)
        assert dep["id"] == 1
        assert post_mock.call_args.args[0].startswith(
            f"{zenodo.ZENODO_SANDBOX_API}/deposit/depositions"
        )


def test_upload_file_puts_to_bucket(tmp_path: Path):
    deposition = {"links": {"bucket": "https://sandbox.zenodo.org/api/files/bucket"}}
    file_path = tmp_path / "file.tar.gz"
    file_path.write_bytes(b"content")

    put_resp = MagicMock()
    put_resp.json.return_value = {"ok": True}
    put_resp.raise_for_status.return_value = None

    with patch("vcfcache.integrations.zenodo.requests.put", return_value=put_resp) as put_mock:
        out = zenodo.upload_file(deposition, file_path, "tok", sandbox=True)
        assert out["ok"] is True
        called_url = put_mock.call_args.args[0]
        assert called_url.endswith("/file.tar.gz")


def test_publish_deposit_posts_publish_link():
    deposition = {"links": {"publish": "https://sandbox.zenodo.org/api/deposit/depositions/1/actions/publish"}}
    post_resp = MagicMock()
    post_resp.json.return_value = {"doi": "10.5072/zenodo.1"}
    post_resp.raise_for_status.return_value = None

    with patch("vcfcache.integrations.zenodo.requests.post", return_value=post_resp) as post_mock:
        out = zenodo.publish_deposit(deposition, "tok", sandbox=True)
        assert out["doi"] == "10.5072/zenodo.1"
        assert post_mock.call_args.args[0] == deposition["links"]["publish"]


def test_resolve_zenodo_alias_queries_keywords_and_honors_sandbox():
    resp = MagicMock()
    resp.ok = True
    resp.json.return_value = {"hits": {"hits": [{"doi": "10.5072/zenodo.12345"}]}}

    with patch("vcfcache.integrations.zenodo.requests.get", return_value=resp) as get_mock:
        doi, alias = zenodo.resolve_zenodo_alias(
            "cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic",
            item_type="caches",
            sandbox=True,
        )
        assert doi == "10.5072/zenodo.12345"
        assert alias == "cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic"

        called_url = get_mock.call_args.args[0]
        called_params = get_mock.call_args.kwargs["params"]
        called_timeout = get_mock.call_args.kwargs["timeout"]
        assert called_url == f"{zenodo.ZENODO_SANDBOX_API}/records/"
        assert "keywords:vcfcache" in called_params["q"]
        assert "keywords:cache" in called_params["q"]
        assert 'keywords:"cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic"' in called_params["q"]
        assert called_timeout == zenodo.DEFAULT_HTTP_TIMEOUT


def test_search_zenodo_records_uses_timeout_and_filters_small_records():
    resp = MagicMock()
    resp.ok = True
    resp.json.return_value = {
        "hits": {
            "hits": [
                {
                    "metadata": {"title": "tiny", "publication_date": "2025-01-01", "keywords": ["vcfcache", "cache"]},
                    "doi": "10.5072/zenodo.tiny",
                    "files": [{"size": 10}],
                },
                {
                    "metadata": {"title": "big", "publication_date": "2025-01-01", "keywords": ["vcfcache", "cache"]},
                    "doi": "10.5072/zenodo.big",
                    "files": [{"size": 2 * 1024 * 1024}],
                },
            ]
        }
    }

    with patch("vcfcache.integrations.zenodo.requests.get", return_value=resp) as get_mock:
        out = zenodo.search_zenodo_records(item_type="caches", sandbox=True, min_size_mb=1.0)
        assert len(out) == 1
        assert out[0]["doi"] == "10.5072/zenodo.big"
        assert get_mock.call_args.kwargs["timeout"] == zenodo.DEFAULT_HTTP_TIMEOUT


def test_resolve_zenodo_alias_raises_on_missing_hits():
    resp = MagicMock()
    resp.ok = True
    resp.json.return_value = {"hits": {"hits": []}}

    with patch("vcfcache.integrations.zenodo.requests.get", return_value=resp):
        with pytest.raises(zenodo.ZenodoError):
            zenodo.resolve_zenodo_alias(
                "cache-hg38-gnomad-4.1joint-AF0100-vep-115.2-basic",
                item_type="caches",
                sandbox=False,
            )


def test_search_zenodo_records_uses_hard_timeout_wrapper(monkeypatch: pytest.MonkeyPatch):
    resp = MagicMock()
    resp.ok = True
    resp.json.return_value = {"hits": {"hits": []}}

    entered: dict[str, bool] = {"ok": False}

    class _NoOpCtx:
        def __enter__(self):
            entered["ok"] = True

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setenv("VCFCACHE_ZENODO_HARD_TIMEOUT", "1")
    with patch("vcfcache.integrations.zenodo._hard_timeout", return_value=_NoOpCtx()):
        with patch("vcfcache.integrations.zenodo.requests.get", return_value=resp):
            zenodo.search_zenodo_records(item_type="caches", sandbox=True, min_size_mb=1.0)

    assert entered["ok"] is True
