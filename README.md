# Sidebot (R Edition)

This is a demonstration of using an LLM to enhance a data dashboard written in Shiny [r-chatbot-cellpanel]().

To run locally, you'll need to create an `.Renviron` file in the repo root with `GEMINI_API_KEY=` followed by a valid GEMINI API key. Or if that environment value is set some other way, you can skip the `.Renviron` file.

Original data was in xlsx format. It was first cleaned up with cell line names obscured, then converted to duckdb with duckdb R package.

## Warnings and limitations

This app sends at least your data schema to a remote LLM. As written, it also permits the LLM to run SQL queries against your data and get the results back. Please keep these facts in mind when dealing with sensitive data. The illustrative example here uses data with cell line names obscured.
